from canonymous.visualization import print_molecule
from canonymous.utils import relabel_molecule
from canonymous.element_properties import ELEMENT_PROPS
from operator import gt, lt, eq
from collections import deque
import networkx as nx


def graph_from_molfile(filename):
    element_symbols, bonds = _parse_molfile(filename)
    element_colors = [ELEMENT_PROPS[s]["element_color"] for s in element_symbols]
    atomic_numbers = [ELEMENT_PROPS[s]["atomic_number"] for s in element_symbols]
    node_labels = range(len(element_symbols))
    graph = nx.Graph()
    graph.add_nodes_from(node_labels)
    graph.add_edges_from(bonds)
    nx.set_node_attributes(graph, dict(zip(node_labels, element_symbols)), "element_symbol")
    nx.set_node_attributes(graph, dict(zip(node_labels, element_colors)), "element_color")
    nx.set_node_attributes(graph, dict(zip(node_labels, atomic_numbers)), "atomic_number")
    nx.set_node_attributes(graph, _find_ring_neighbors(graph), "ring_neighbors")
    nx.set_node_attributes(graph, 0, "partition")
    _add_fingerprint_attribute(graph)
    return graph

def _parse_molfile(filename):
    with open(filename) as f:
        lines = [l.rstrip().split(" ") for l in f]
    lines = [[value for value in line if value != ""] for line in lines]
    atom_count = int(lines[5][3])
    bond_count = int(lines[5][4])
    atom_block_offset = 7
    bond_block_offset = atom_block_offset + atom_count + 2
    element_symbols = [l[3] for l in lines[atom_block_offset:atom_block_offset + atom_count]]
    bonds = [(int(l[4]) - 1, int(l[5]) - 1)
             for l in lines[bond_block_offset:bond_block_offset + bond_count]]    # make bond-indices zero-based
    return element_symbols, bonds

def _add_fingerprint_attribute(m):
    """Assign a fingerprint of the following format to each atom:
    <atomic number><indices ring neighbors sorted in non-ascending order>"""
    atomic_numbers = list(nx.get_node_attributes(m, "atomic_number").values())
    ring_neighbors = list(nx.get_node_attributes(m, "ring_neighbors").values())
    ring_neighbors = [sorted((list(rn)), reverse=True) for rn in ring_neighbors]
    fingerprints = [str(a) + "".join(map(str, r)) for a, r in zip(atomic_numbers, ring_neighbors)]
    nx.set_node_attributes(m, dict(zip(range(m.number_of_nodes()), fingerprints)), "fingerprint")

def sort_molecule_by_attribute(m, attribute):
    '''Sort atoms lexicographically by attribute.'''
    attr_sequence = [_attribute_sequence(atom, m, attribute) for atom in m]
    idcs = list(range(m.number_of_nodes()))
    attr_with_idcs = [(i, j) for i, j in zip(attr_sequence, idcs)] # [(A, 0), (C, 1), (B, 2)]
    sorted_attr, idcs_sorted_by_attr = zip(*sorted(attr_with_idcs)) # (A, B, C), (0, 2, 1)
    return relabel_molecule(m, idcs_sorted_by_attr, idcs)

def _attribute_sequence(atom, m, attribute):
    attr_atom = m.nodes[atom][attribute]
    attr_neighbors = sorted([m.nodes[n][attribute] for n in m.neighbors(atom)], reverse=True)
    return [attr_atom] + attr_neighbors

def partition_molecule_by_attribute(m, attribute):
    current_partition = 0
    for i in range(m.number_of_nodes() - 1):
        j = i + 1
        atomic_numbers_i = _attribute_sequence(i, m, attribute)
        atomic_numbers_j = _attribute_sequence(j, m, attribute)
        if (atomic_numbers_i != atomic_numbers_j): current_partition += 1
        m.nodes[j]["partition"] = current_partition
    return m

def partition_molecule_recursively(m, show_steps=False):
    m_sorted = sort_molecule_by_attribute(m, "partition")
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
        if (partitions_i != partitions_j):
            current_partition += 1
        updated_partitions.append(current_partition)
    if current_partitions == updated_partitions:
        return m_sorted
    nx.set_node_attributes(m_sorted,
                           dict(zip(range(n_nodes), updated_partitions)),
                           "partition")
    return partition_molecule_recursively(m_sorted, show_steps=show_steps)

def assign_canonical_labels(m, root_idx, traversal_priorities=[lt, gt, eq], show_traversal_order=False):
    partitions = m.nodes.data("partition")
    lut = _create_partition_lut(m)
    atom_queue = [root_idx]
    canonical_idcs = {}
    nx.set_node_attributes(m, False, "explored")

    while atom_queue:
        a = atom_queue.pop()
        if m.nodes[a]["explored"]:
            continue
        a_canon = lut[partitions[a]].pop()
        canonical_idcs[a] = a_canon
        if show_traversal_order:
            print(f"Current atom index: {a}.\tRe-labeling to {a_canon}.")
        neighbors = list(m.neighbors(a))
        neighbor_traversal_order = []
        for priority in traversal_priorities:
            neighbors_this_priority = [n for n in neighbors
                                       if priority(partitions[a], partitions[n])]
            neighbor_traversal_order.extend(sorted(neighbors_this_priority))

        m.nodes[a]["explored"] = True
        for n in neighbor_traversal_order:
            atom_queue.insert(0, n)

    nx.set_node_attributes(m, False, "explored")
    return canonical_idcs

def _create_partition_lut(m):
    """Look-up-table of atom indices in partitions."""
    partitions = set(sorted([v for k, v in m.nodes.data("partition")]))
    partition_lut = {p:set() for p in partitions}
    for a in m:
        partition_lut[m.nodes[a]["partition"]].add(a)
    partition_lut.update((k, sorted(list(v), reverse=True)) for k, v in partition_lut.items())
    return partition_lut

def canonicalize_molecule(m, root_idx=0):
    m_sorted_by_fingerprint = sort_molecule_by_attribute(m, "fingerprint")
    m_partitioned_by_fingerprint = partition_molecule_by_attribute(m_sorted_by_fingerprint, "fingerprint")
    m_partitioned = partition_molecule_recursively(m_partitioned_by_fingerprint, show_steps=False)
    canonical_idcs = assign_canonical_labels(m_partitioned, root_idx)
    return nx.relabel_nodes(m_partitioned, canonical_idcs, copy=True)

def _find_atomic_cycles(m):
    """Find atomic cycles in a graph.

    An atomic cycle is a generalization of a chordless cycle, such that the
    chord can be longer than one edge.

    Adaptation of the `find_large_atomic_cycle` algorithm from figure 6 in
    DOI: 10.1080/09540091.2012.664122. In constrast to the original
    implementation, the present implementation returns *all* cycles and uses a
    fixed lambda threshold of size 3.
    """
    root_node = 0
    outer_node_queue = deque([root_node])    # Q
    outer_visited_nodes = set([root_node])    # S
    visited_edges = set()    # T
    atomic_cycles = []

    while outer_node_queue:    # outer BFS
        outer_node = outer_node_queue.popleft()    # a
        for outer_neighbor in m.neighbors(outer_node):    # b
            if outer_neighbor not in outer_visited_nodes:
                outer_node_queue.append(outer_neighbor)
                outer_visited_nodes.add(outer_neighbor)
                visited_edges.add((outer_node, outer_neighbor))
                visited_edges.add((outer_neighbor, outer_node))
                continue

            # A cycle has been detected.
            inner_node_queue = deque([outer_neighbor])    # I
            inner_visited_nodes = set([outer_neighbor])    # U
            parent_nodes = {outer_neighbor: -1}    # P

            while inner_node_queue:    # inner BFS
                inner_node = inner_node_queue.popleft()    # c
                for inner_neighbor in m.neighbors(inner_node):    # d
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
                    cycle = set()    # Y
                    while parent_nodes[inner_neighbor] != -1:
                        cycle_node, inner_neighbor = inner_neighbor, parent_nodes[inner_neighbor]
                        cycle.add(cycle_node)
                    cycle.add(inner_neighbor)
                    if len(cycle) >= 3:
                        atomic_cycles.append(cycle)
                    inner_node_queue.clear()
                    break
            visited_edges.add((outer_node, outer_neighbor))
            visited_edges.add((outer_neighbor, outer_node))
    return atomic_cycles

def _find_ring_neighbors(m):
    ring_neighbors = {atom: {-1} for atom in m.nodes}
    for ring in _find_atomic_cycles(m):
        for atom in ring:
            ring_neighbors[atom].update(ring)
    return ring_neighbors

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
