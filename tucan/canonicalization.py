from tucan.graph_utils import sort_molecule_by_attribute, attribute_sequence
import networkx as nx
from igraph import Graph as iGraph
from itertools import pairwise


def partition_molecule_by_attribute(m, attribute):
    m_sorted = sort_molecule_by_attribute(m, attribute)
    sorted_indices = sorted(m_sorted)
    attribute_sequences = [
        attribute_sequence(a, m_sorted, attribute) for a in sorted_indices
    ]
    updated_partitions = [0]
    for i, j in pairwise(attribute_sequences):
        current_partition = updated_partitions[-1]
        if i != j:
            current_partition += 1
        updated_partitions.append(current_partition)
    nx.set_node_attributes(
        m_sorted, dict(zip(sorted_indices, updated_partitions)), "partition"
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


def assign_canonical_labels(m):
    """Canonicalize node-labels of a graph.

    The canonical labels are computed using the igraph [1] implementation of
    the "bliss" algorithm [2].

    Returns
    -------
    dict
        From old labels (keys) to canonical labels (values).

    References
    ----------
    [1] https://igraph.org
    [2] https://doi.org/10.1137/1.9781611972870.13
    """

    m_igraph = iGraph.from_networkx(m)
    old_labels = m_igraph.vs["_nx_name"]
    partitions = m_igraph.vs["partition"]
    canonical_labels = m_igraph.canonical_permutation(color=partitions)

    return dict(zip(old_labels, canonical_labels))


def _add_invariant_code(m, invariant_code_definitions):
    """Assign an invariant code to each atom (mutates graph)."""
    invariant_codes = [
        tuple(
            attributes[icd["key"]]
            if (default_value := icd["default_value"]) is None
            else attributes.get(icd["key"], default_value)
            for icd in invariant_code_definitions
        )
        for _, attributes in m.nodes(data=True)
    ]

    nx.set_node_attributes(
        m, dict(zip(list(m.nodes), invariant_codes)), "invariant_code"
    )


def _invariant_code_definition(attribute_key, default_value=None):
    """
    Returns an invariant code definition ("icd") to be used in the _add_invariant_code function.
    Parameters
    ----------
    attribute_key Node attribute key
    default_value Default value to be used if the node attribute does not exist. Use None to indicate that the attribute is mandatory.
    """
    return {"key": attribute_key, "default_value": default_value}


def canonicalize_molecule(m):
    invariant_code_definitions = [
        _invariant_code_definition("atomic_number"),
        _invariant_code_definition("mass", 0),
        _invariant_code_definition("rad", 0),
    ]
    _add_invariant_code(m, invariant_code_definitions)

    m_partitioned_by_invariant_code = partition_molecule_by_attribute(
        m, "invariant_code"
    )
    m_refined = list(refine_partitions(m_partitioned_by_invariant_code))
    m_partitioned = m_refined[-1] if m_refined else m_partitioned_by_invariant_code
    canonical_labels = assign_canonical_labels(m_partitioned)
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
