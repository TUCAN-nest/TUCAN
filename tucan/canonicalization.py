from operator import gt, lt, eq
from tucan.graph_utils import sort_molecule_by_attribute, attribute_sequence
import networkx as nx


def partition_molecule_by_attribute(m, attribute, include_neighbors=True):
    m_sorted = sort_molecule_by_attribute(m, attribute)
    updated_partitions = [0]
    n_nodes = m.number_of_nodes()
    for i in range(n_nodes - 1):
        j = i + 1
        attributes_i = (
            attribute_sequence(i, m_sorted, attribute)
            if include_neighbors
            else m_sorted.nodes[i][attribute]
        )
        attributes_j = (
            attribute_sequence(j, m_sorted, attribute)
            if include_neighbors
            else m_sorted.nodes[j][attribute]
        )
        current_partition = updated_partitions[-1]
        if attributes_i != attributes_j:
            current_partition += 1
        updated_partitions.append(current_partition)
    nx.set_node_attributes(
        m_sorted, dict(zip(range(n_nodes), updated_partitions)), "partition"
    )
    return m_sorted


def refine_partitions(m):
    current_partitions = list(
        nx.get_node_attributes(
            sort_molecule_by_attribute(m, "partition"), "partition"
        ).values()
    )
    m_refined = partition_molecule_by_attribute(m, "partition")
    refined_partitions = list(nx.get_node_attributes(m_refined, "partition").values())

    while current_partitions != refined_partitions:
        yield m_refined
        current_partitions = refined_partitions
        m_refined = partition_molecule_by_attribute(m_refined, "partition")
        refined_partitions = list(
            nx.get_node_attributes(m_refined, "partition").values()
        )


def assign_canonical_labels(
    m,
    root_idx,
    traversal_priorities=[lt, gt, eq],
    enforce_deterministic_branching=False,
    show_traversal_order=False,
):
    partitions = m.nodes.data("partition")
    labels_by_partition = _labels_by_partition(m)
    atom_queue = [root_idx]
    canonical_labels = {}
    nx.set_node_attributes(m, False, "explored")

    if enforce_deterministic_branching:
        partitions_to_permute = set()

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
        for priority in reversed(
            traversal_priorities
        ):  # reverse to preserve order of traversal priorities in queue
            neighbors_this_priority = [
                n for n in neighbors if priority(partitions[a], partitions[n])
            ]
            neighbor_traversal_order.extend(sorted(neighbors_this_priority))

        m.nodes[a]["explored"] = True

        if enforce_deterministic_branching:
            neighbor_partitions = [
                partitions[a]
                for a in neighbor_traversal_order
                if (len(list(m.neighbors(a))) > 1) and (not m.nodes[a]["explored"])
            ]  # non-terminal neighbors only
            duplicate_partitions = [
                p for p in neighbor_partitions if neighbor_partitions.count(p) > 1
            ]
            partitions_to_permute.update(duplicate_partitions)

        for n in neighbor_traversal_order:
            atom_queue.insert(0, n)

    nx.set_node_attributes(m, False, "explored")

    if enforce_deterministic_branching:
        n_partitions = len(set(dict(partitions).values()))
        n_partitions_to_permute = len(partitions_to_permute)
        assert (
            not partitions_to_permute
        ), f"Couldn't branch determinstically to {n_partitions_to_permute}/{n_partitions} partition(s): {partitions_to_permute}."

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


def _add_invariant_code(m):
    """Assign an invariant code to each atom (mutates graph)."""
    atomic_numbers = list(nx.get_node_attributes(m, "atomic_number").values())
    invariant_codes = list(map(str, atomic_numbers))
    nx.set_node_attributes(
        m, dict(zip(range(m.number_of_nodes()), invariant_codes)), "invariant_code"
    )

def canonicalize_molecule(m, root_idx=0):
    _add_invariant_code(m)
    m_partitioned_by_invariant_code = partition_molecule_by_attribute(
        m, "invariant_code"
    )
    m_refined = list(refine_partitions(m_partitioned_by_invariant_code))
    m_partitioned = m_refined[-1] if m_refined else m_partitioned_by_invariant_code
    canonical_labels = assign_canonical_labels(m_partitioned, root_idx)
    return nx.relabel_nodes(m_partitioned, canonical_labels, copy=True)


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
