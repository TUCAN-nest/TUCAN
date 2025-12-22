from tucan.graph_attributes import PARTITION, INVARIANT_CODE
from tucan.graph_utils import get_attribute_sequences
import networkx as nx
from typing import Generator, Any
from collections import Counter


def get_partitions_from_attribute(
    attributes: dict[int, Any], neighbors: dict[int, tuple[int]]
) -> dict[int, int]:
    # Node degree (i.e., number of neighbors) is encoded in length of individual attribute sequences.
    attr_seqs = get_attribute_sequences(attributes, neighbors)
    unique_attr_seqs = sorted(set(attr_seqs))
    unique_attr_seqs_to_partitions = {
        attr_seq: partition for partition, attr_seq in enumerate(unique_attr_seqs)
    }

    return dict(
        zip(
            attributes.keys(),
            [unique_attr_seqs_to_partitions[attr_seq] for attr_seq in attr_seqs],
            strict=True,
        )
    )


def set_partitions(
    m: nx.Graph, partitions: dict[int, int], copy: bool = True
) -> nx.Graph:
    m_partitioned = m.copy() if copy else m

    for atom, partition in partitions.items():
        m_partitioned.nodes[atom][
            PARTITION
        ] = partition  # Setting partition attribute directly is faster than using nx.set_node_attributes() for batch update.

    m_partitioned.graph["n_partitions"] = len(set(partitions.values()))

    return m_partitioned


def partition_molecule_by_attribute(m: nx.Graph, attribute: str) -> nx.Graph:
    attributes = dict(nx.get_node_attributes(m, attribute))
    neighbors = {node: tuple(m[node]) for node in m}
    partitions = get_partitions_from_attribute(attributes, neighbors)

    return set_partitions(m, partitions)


def _refine_partitions(m: nx.Graph) -> Generator[dict[int, int], None, None]:
    partitions = dict(nx.get_node_attributes(m, PARTITION))
    neighbors = {node: tuple(m[node]) for node in m}

    while True:
        yield partitions
        refined = get_partitions_from_attribute(partitions, neighbors)
        if refined == partitions:
            return
        partitions = refined


def refine_partitions(m: nx.Graph) -> Generator[nx.Graph, None, None]:
    return (set_partitions(m, partitions) for partitions in _refine_partitions(m))


def get_target_partition(m: nx.Graph) -> int:
    partitions = nx.get_node_attributes(m, PARTITION).values()
    partition_sizes = Counter(sorted(partitions))

    return max(partition_sizes, key=partition_sizes.get)  # type: ignore


def get_refinement_tree_node_children(m: nx.Graph) -> Generator[nx.Graph, None, None]:
    n_partitions = m.graph["n_partitions"]
    target_partition = get_target_partition(m)

    for atom, partition in m.nodes(data=PARTITION):  # type: ignore
        if partition != target_partition:
            continue

        # Split target partition.
        m_artificially_split = (
            nx.Graph()
        )  # Instantiate a new graph instead of using m.copy() to improve performance.
        m_artificially_split.add_nodes_from(m.nodes(data=True))
        m_artificially_split.add_edges_from(m.edges(data=True))
        m_artificially_split.nodes[atom][
            PARTITION
        ] = n_partitions  # Partitions are zero-based.
        m_artificially_split.graph["n_partitions"] = n_partitions + 1

        m_artificially_refined = set_partitions(
            m_artificially_split,
            list(_refine_partitions(m_artificially_split))[-1],
            False,
        )

        yield m_artificially_refined


def filter_out_automorphisms(ms: list[nx.Graph]) -> list[nx.Graph]:
    # Caution: Mutates `ms` in-place.
    filtered_ms = set()
    labelings = set()

    for m in ms:
        partitions = list(nx.get_node_attributes(m, PARTITION).values())
        _labeling = []
        for edge in m.edges():
            edge_partitions = [partitions[n] for n in edge]
            edge_partitions.sort()
            _labeling.append(tuple(edge_partitions))

        _labeling.sort()
        labeling = tuple(_labeling)

        if labeling in labelings:
            continue

        labelings.add(labeling)
        filtered_ms.add(m)

    return list(filtered_ms)


def partitioning_is_discrete(m):
    return m.graph["n_partitions"] == m.number_of_nodes()


def get_refinement_tree_levels(
    m: nx.Graph, filter_automorphisms: bool = True
) -> Generator[list[nx.Graph], None, None]:
    """
    Build BFS refinement-tree and yield each level.
    """
    parents = [m]

    while parents:
        yield parents
        if all(map(partitioning_is_discrete, parents)):
            return

        children = [
            child
            for parent in parents
            for child in get_refinement_tree_node_children(parent)
        ]
        parents = (
            filter_out_automorphisms(children) if filter_automorphisms else children
        )


def get_canonical_molecule(ms: list[nx.Graph]) -> nx.Graph:
    m_canonical = ms[0]
    canonical_labeling = [[0, 0]]

    for m in ms:
        m_relabeled_by_partition = nx.relabel_nodes(
            m,
            dict(m.nodes(data=PARTITION)),
            copy=True,
        )
        labeling = sorted([sorted(edge) for edge in m_relabeled_by_partition.edges()])

        if labeling > canonical_labeling:
            m_canonical = m_relabeled_by_partition
            canonical_labeling = labeling

    return m_canonical


def canonicalize_molecule(m: nx.Graph) -> nx.Graph:
    m_partitioned_by_invariant_code = partition_molecule_by_attribute(m, INVARIANT_CODE)
    m_refined = list(refine_partitions(m_partitioned_by_invariant_code))[-1]
    ms_discrete_partitionings = list(get_refinement_tree_levels(m_refined))[-1]

    return get_canonical_molecule(ms_discrete_partitionings)
