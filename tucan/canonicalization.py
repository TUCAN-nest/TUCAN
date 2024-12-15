from tucan.graph_attributes import PARTITION, INVARIANT_CODE
from tucan.graph_utils import attribute_sequence
import networkx as nx
from typing import Generator
from collections import Counter


def partition_molecule_by_attribute(
    m: nx.Graph, attribute: str, copy: bool = True
) -> nx.Graph:
    # Node degree (i.e., number of neighbors) is encoded in length of individual attribute sequences.
    attr_seqs = [attribute_sequence(m, atom, attribute) for atom in m]
    unique_attr_seqs = sorted(set(attr_seqs))
    unique_attr_seqs_to_partitions = dict(
        zip(unique_attr_seqs, range(len(unique_attr_seqs)))
    )
    partitions = [unique_attr_seqs_to_partitions[attr_seq] for attr_seq in attr_seqs]

    m_partitioned = m.copy() if copy else m
    nx.set_node_attributes(
        m_partitioned, dict(zip(list(m_partitioned), partitions)), PARTITION
    )
    m_partitioned.graph["n_partitions"] = len(unique_attr_seqs)

    return m_partitioned


def partitioning_is_discrete(m):
    return m.graph["n_partitions"] == m.number_of_nodes()


def refine_partitions(m: nx.Graph) -> Generator[nx.Graph, None, None]:
    n_partitions = m.graph["n_partitions"]
    m_refined = partition_molecule_by_attribute(m, PARTITION, copy=False)
    n_partitions_refined = m_refined.graph["n_partitions"]

    if n_partitions == n_partitions_refined:
        # No more refinement possible.
        yield m_refined
        return

    yield from refine_partitions(m_refined)


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
        m_artificially_split = m.copy()
        nx.set_node_attributes(
            m_artificially_split,
            {atom: n_partitions},  # Partitions are zero-based.
            PARTITION,
        )
        m_artificially_refined = list(refine_partitions(m_artificially_split))[-1]

        yield m_artificially_refined


def filter_out_automorphisms(ms: list[nx.Graph]) -> list[nx.Graph]:
    # Caution: Mutates `ms` in-place.
    filtered_ms = set()
    labelings = set()

    for m in ms:
        m_relabeled_by_partition = nx.relabel_nodes(
            m,
            dict(zip(list(m), nx.get_node_attributes(m, PARTITION).values())),
        )
        labeling = tuple(
            sorted([tuple(sorted(edge)) for edge in m_relabeled_by_partition.edges()])
        )

        if labeling in labelings:
            continue

        labelings.add(labeling)
        filtered_ms.add(m)

    return list(filtered_ms)


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
            dict(zip(list(m), nx.get_node_attributes(m, PARTITION).values())),
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
