from canonymous.visualization import print_molecule
from canonymous.element_properties import ELEMENT_PROPS
from operator import gt, lt, eq
from collections import deque, Counter
import networkx as nx
import random


def graph_from_molfile(filename):
    element_symbols, bonds = _parse_molfile(filename)
    element_colors = [ELEMENT_PROPS[s]["element_color"] for s in element_symbols]
    atomic_numbers = [ELEMENT_PROPS[s]["atomic_number"] for s in element_symbols]
    node_labels = range(len(element_symbols))
    graph = nx.Graph()
    graph.add_nodes_from(node_labels)
    graph.add_edges_from(bonds)
    nx.set_node_attributes(
        graph, dict(zip(node_labels, element_symbols)), "element_symbol"
    )
    nx.set_node_attributes(
        graph, dict(zip(node_labels, element_colors)), "element_color"
    )
    nx.set_node_attributes(
        graph, dict(zip(node_labels, atomic_numbers)), "atomic_number"
    )
    nx.set_node_attributes(graph, _cycle_memberships(graph), "cycle_membership")
    nx.set_node_attributes(graph, 0, "partition")
    _add_invariant_code(graph)
    return graph


def _parse_molfile(filename):
    with open(filename) as f:
        lines = [l.rstrip().split(" ") for l in f]
    lines = [[value for value in line if value != ""] for line in lines]
    atom_count = int(lines[5][3])
    bond_count = int(lines[5][4])
    atom_block_offset = 7
    bond_block_offset = atom_block_offset + atom_count + 2
    element_symbols = [
        l[3] for l in lines[atom_block_offset : atom_block_offset + atom_count]
    ]
    bonds = [
        (int(l[4]) - 1, int(l[5]) - 1)
        for l in lines[bond_block_offset : bond_block_offset + bond_count]
    ]  # make bond-indices zero-based
    return element_symbols, bonds


def _cycle_memberships(m):
    cycle_memberships = {node: {0} for node in m.nodes}
    for cycle_id, cycle in enumerate(_find_atomic_cycles(m)):
        for node in cycle:
            cycle_memberships[node].add(cycle_id + 1)
    return cycle_memberships


def _find_atomic_cycles(m):
    """Find atomic cycles in a graph.

    An atomic cycle is a generalization of a chordless cycle, such that the
    chord can be longer than one edge. Adaptation of the
    `find_large_atomic_cycle` algorithm from Figure 6 in [1]. In contrast to
    [1] the present implementation returns all cycles and uses a fixed lambda
    threshold of size 3.

    References
    ----------
    [1] DOI: 10.1080/09540091.2012.664122
    """
    root_node = 0
    outer_node_queue = deque([root_node])  # Q
    outer_visited_nodes = set([root_node])  # S
    visited_edges = set()  # T
    atomic_cycles = []

    while outer_node_queue:  # outer BFS
        outer_node = outer_node_queue.popleft()  # a
        for outer_neighbor in m.neighbors(outer_node):  # b
            if outer_neighbor not in outer_visited_nodes:
                outer_node_queue.append(outer_neighbor)
                outer_visited_nodes.add(outer_neighbor)
                visited_edges.add((outer_node, outer_neighbor))
                visited_edges.add((outer_neighbor, outer_node))
                continue

            # A cycle has been detected.
            inner_node_queue = deque([outer_neighbor])  # I
            inner_visited_nodes = set([outer_neighbor])  # U
            parent_nodes = {outer_neighbor: -1}  # P

            while inner_node_queue:  # inner BFS
                inner_node = inner_node_queue.popleft()  # c
                for inner_neighbor in m.neighbors(inner_node):  # d
                    if inner_neighbor in inner_visited_nodes:
                        continue
                    if (inner_neighbor, inner_node) not in visited_edges:
                        continue
                    parent_nodes[inner_neighbor] = inner_node
                    inner_node_queue.append(inner_neighbor)
                    inner_visited_nodes.add(inner_neighbor)

                    if inner_neighbor != outer_node:
                        continue
                    # An atomic cycle has been detected.
                    cycle = set()  # Y
                    while parent_nodes[inner_neighbor] != -1:
                        cycle_node, inner_neighbor = (
                            inner_neighbor,
                            parent_nodes[inner_neighbor],
                        )
                        cycle.add(cycle_node)
                    cycle.add(inner_neighbor)
                    if len(cycle) >= 3:
                        atomic_cycles.append(cycle)
                    inner_node_queue.clear()
                    break
            visited_edges.add((outer_node, outer_neighbor))
            visited_edges.add((outer_neighbor, outer_node))
    return atomic_cycles


def _add_invariant_code(m):
    """Assign an invariant code to each atom.

    Invariant code is formatted as follows:
    <atomic number>-<comma-separated list of fundamental cycles>
    """
    atomic_numbers = list(nx.get_node_attributes(m, "atomic_number").values())
    cycle_memberships = list(nx.get_node_attributes(m, "cycle_membership").values())
    cycle_memberships = [sorted((list(cm)), reverse=True) for cm in cycle_memberships]
    invariant_codes = [
        f"{str(a)}-{','.join(map(str, c))}"
        for a, c in zip(atomic_numbers, cycle_memberships)
    ]
    nx.set_node_attributes(
        m, dict(zip(range(m.number_of_nodes()), invariant_codes)), "invariant_code"
    )


def partition_molecule_by_attribute(m, attribute, include_neighbors=True):
    m_sorted = _sort_molecule_by_attribute(m, attribute)
    current_partition = 0
    for i in range(m_sorted.number_of_nodes() - 1):
        j = i + 1
        attributes_i = (
            _attribute_sequence(i, m_sorted, attribute)
            if include_neighbors
            else m_sorted.nodes[i][attribute]
        )
        attributes_j = (
            _attribute_sequence(j, m_sorted, attribute)
            if include_neighbors
            else m_sorted.nodes[j][attribute]
        )
        if attributes_i != attributes_j:
            current_partition += 1
        m_sorted.nodes[j]["partition"] = current_partition
    return m_sorted


def _sort_molecule_by_attribute(m, attribute):
    """Sort atoms by attribute."""
    attr_sequence = [_attribute_sequence(atom, m, attribute) for atom in m]
    labels = list(range(m.number_of_nodes()))
    attr_with_labels = [
        (i, j) for i, j in zip(attr_sequence, labels)
    ]  # [(A, 0), (C, 1), (B, 2)]
    sorted_attr, labels_sorted_by_attr = zip(
        *sorted(attr_with_labels)
    )  # (A, B, C), (0, 2, 1)
    return _relabel_molecule(m, labels_sorted_by_attr, labels)


def _attribute_sequence(atom, m, attribute):
    attr_atom = m.nodes[atom][attribute]
    attr_neighbors = sorted(
        [m.nodes[n][attribute] for n in m.neighbors(atom)], reverse=True
    )
    return [attr_atom] + attr_neighbors


def _relabel_molecule(m, old_labels, new_labels):
    """Relabel the atoms of a molecular graph."""
    m_relabeled = nx.relabel_nodes(m, dict(zip(old_labels, new_labels)))
    # In the NetworkX Graph datastructure, the relabeled nodes don't occur in
    # increasing order yet. This is why we change the node order now.
    m_sorted = nx.Graph()
    m_sorted.add_nodes_from(sorted(m_relabeled.nodes(data=True)))
    m_sorted.add_edges_from(m_relabeled.edges(data=True))
    return m_sorted


def partition_molecule_recursively(m, show_steps=False):
    m_sorted = _sort_molecule_by_attribute(m, "partition")
    if show_steps:
        print_molecule(m_sorted, "refined partitions")
    current_partitions = list(nx.get_node_attributes(m, "partition").values())
    updated_partitions = [0]
    n_nodes = m.number_of_nodes()
    for i in range(n_nodes - 1):
        j = i + 1
        partitions_i = _attribute_sequence(i, m_sorted, "partition")
        partitions_j = _attribute_sequence(j, m_sorted, "partition")
        current_partition = updated_partitions[-1]
        if partitions_i != partitions_j:
            current_partition += 1
        updated_partitions.append(current_partition)
    if current_partitions == updated_partitions:
        return m_sorted
    nx.set_node_attributes(
        m_sorted, dict(zip(range(n_nodes), updated_partitions)), "partition"
    )
    return partition_molecule_recursively(m_sorted, show_steps=show_steps)


def assign_canonical_labels(
    m, root_idx, traversal_priorities=[lt, gt, eq], show_traversal_order=False
):
    partitions = m.nodes.data("partition")
    labels_by_partition = _labels_by_partition(m)
    atom_queue = [root_idx]
    canonical_labels = {}
    nx.set_node_attributes(m, False, "explored")

    while atom_queue:
        a = atom_queue.pop()
        if m.nodes[a]["explored"]:
            continue
        a_canon = labels_by_partition[partitions[a]].pop()
        canonical_labels[a] = a_canon
        if show_traversal_order:
            print(f"Current atom index: {a}.\tRe-labeling to {a_canon}.")
        neighbors = list(m.neighbors(a))
        neighbor_traversal_order = []
        for priority in traversal_priorities:
            neighbors_this_priority = [
                n for n in neighbors if priority(partitions[a], partitions[n])
            ]
            neighbor_traversal_order.extend(sorted(neighbors_this_priority))

        m.nodes[a]["explored"] = True
        for n in neighbor_traversal_order:
            atom_queue.insert(0, n)

    nx.set_node_attributes(m, False, "explored")
    return canonical_labels


def _labels_by_partition(m):
    """Create dictionary of partitions to atom labels."""
    partitions = set(sorted([v for k, v in m.nodes.data("partition")]))
    labels_by_partition = {p: set() for p in partitions}
    for a in m:
        labels_by_partition[m.nodes[a]["partition"]].add(a)
    labels_by_partition.update(
        (k, sorted(list(v), reverse=True)) for k, v in labels_by_partition.items()
    )
    return labels_by_partition


def canonicalize_molecule(m, root_idx=0):
    m_partitioned_by_invariant_code = partition_molecule_by_attribute(
        m, "invariant_code"
    )
    m_partitioned = partition_molecule_recursively(
        m_partitioned_by_invariant_code, show_steps=False
    )
    canonical_labels = assign_canonical_labels(m_partitioned, root_idx)
    return nx.relabel_nodes(m_partitioned, canonical_labels, copy=True)


def serialize_molecule(m):
    """Serialize a molecule."""
    serialization = _sum_formula(m)
    for edge in sorted([sorted(edge) for edge in m.edges()]):
        serialization += f"/{edge[0]}-{edge[1]}"
    return serialization


def _sum_formula(m):
    """Write sum formula of a molecule.

    Elements occur in the following order: C, H, other elements in alphabetic order.
    """
    element_counts = Counter(nx.get_node_attributes(m, "element_symbol").values())
    element_counts = {
        k: (v if v > 1 else "") for k, v in element_counts.items()
    }  # remove counts of 1 since those are implicit in sum formula
    sum_formula = ""
    for element in ["C", "H"]:
        count = element_counts.pop(element, None)
        if count:
            sum_formula += f"{element}{count}"
    for k, v in dict(sorted(element_counts.items())).items():
        sum_formula += f"{k}{v}"
    return sum_formula


def permute_molecule(m, random_seed=1.0):
    """Randomly permute the atom-labels of a molecular graph.

    Parameters
    ----------
    random_seed: float
        In [0.0, 1.0).
    """
    labels = m.nodes()
    permuted_labels = list(range(m.number_of_nodes()))
    # Enforce permutation for graphs with at least 2 edges that aren't fully connected (i.e., complete).
    enforce_permutation = m.number_of_edges() > 1 and nx.density(m) != 1
    random.seed(
        random_seed
    )  # subsequent calls of random.shuffle(x[, random]) will now use fixed sequence of values for `random` parameter
    random.shuffle(permuted_labels)
    m_permu = _relabel_molecule(m, permuted_labels, labels)
    if enforce_permutation:
        while m.edges == m_permu.edges:
            random.shuffle(permuted_labels)
            m_permu = _relabel_molecule(m, permuted_labels, labels)
    return m_permu


# def bfs_molecule(m, root_idx):
#     """Breadth-first search over atoms.
#     Note that NetworkX provides the same algorithm in `dfs_edges()`.
#     This (re-)implementation allows for controlling the branching behavior
#     during the molecule traversal.
#     m: NetworkX graph.
#     root_idx: atom at which to start traversal.
#     """
#     m.nodes[root_idx]["explored"] = True
#     atom_queue = deque([root_idx])
#     while atom_queue:
#         a = atom_queue.popleft()
#         for n in m.neighbors(a):
#             if m.nodes[n]["explored"]:
#                 continue
#             yield (a, n)
#             m.nodes[n]["explored"] = True
#             atom_queue.append(n)

# def dfs_molecule(m, root_idx):
#     """Depth-first search over atoms.
#     Note that NetworkX provides the same algorithm in `bfs_edges()`.
#     This (re-)implementation allows for controlling the branching behavior
#     during the molecule traversal.
#     m: NetworkX graph.
#     root_idx: atom at which to start traversal.
#     """
#     m.nodes[root_idx]["explored"] = True
#     for n_idx in m.neighbors(root_idx):
#         if m.nodes[n_idx]["explored"]:
#             continue
#         yield (root_idx, n_idx)
#         yield from dfs_molecule(m, n_idx)

# def edge_dfs_molecule(m, root_idx):
#     """Depth-first search over edges.
#     Note that NetworkX provides the same algorithm in `edge_dfs ()`.
#     This (re-)implementation allows for controlling the branching behavior
#     during the molecule traversal.
#     m: NetworkX graph.
#     root_idx: atom at which to start traversal.
#     """
#     visited_edges = set()
#     visited_nodes = set()
#     edges = {}

#     nodes = list(m.nbunch_iter(root_idx))
#     for start_node in nodes:
#         stack = [start_node]
#         while stack:
#             current_node = stack[-1]
#             if current_node not in visited_nodes:
#                 edges[current_node] = iter(m.edges(current_node))
#                 visited_nodes.add(current_node)

#             try:
#                 edge = next(edges[current_node])
#             except StopIteration:
#                 # No more edges from the current node.
#                 stack.pop()
#             else:
#                 edgeid = (frozenset(edge[:2]),) + edge[2:]
#                 if edgeid not in visited_edges:
#                     visited_edges.add(edgeid)
#                     # Mark the traversed "to" node as to-be-explored.
#                     stack.append(edge[1])
#                     yield edge
