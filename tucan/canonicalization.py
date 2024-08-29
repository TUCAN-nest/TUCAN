from tucan.graph_attributes import PARTITION, INVARIANT_CODE
from tucan.graph_utils import attribute_sequence
import networkx as nx
from typing import Generator
from collections import Counter, deque


def partition_molecule_by_attribute(m: nx.Graph, attribute: str) -> nx.Graph:
    # Node degree (i.e., number of neighbors) is encoded in length of individual attribute sequences.
    attr_seqs = [attribute_sequence(m, atom, attribute) for atom in m]
    unique_attr_seqs = sorted(set(attr_seqs))
    unique_attr_seqs_to_partitions = dict(
        zip(unique_attr_seqs, range(len(unique_attr_seqs)))
    )
    partitions = [unique_attr_seqs_to_partitions[attr_seq] for attr_seq in attr_seqs]

    m_partitioned = m.copy()
    nx.set_node_attributes(
        m_partitioned, dict(zip(list(m_partitioned), partitions)), PARTITION
    )

    return m_partitioned


def get_number_of_partitions(m: nx.Graph) -> int:
    return len(set(nx.get_node_attributes(m, PARTITION).values()))


def partitioning_is_discrete(m):
    return get_number_of_partitions(m) == m.number_of_nodes()


def refine_partitions(m: nx.Graph) -> Generator[nx.Graph, None, None]:
    m_refined = partition_molecule_by_attribute(m, PARTITION)

    if get_number_of_partitions(m_refined) == get_number_of_partitions(m):
        # No more refinement possible.
        yield m_refined
        return

    yield from refine_partitions(m_refined)


def get_refinement_tree_node_children(m: nx.Graph) -> Generator[nx.Graph, None, None]:
    n_partitions = get_number_of_partitions(m)
    partitions = nx.get_node_attributes(m, PARTITION)
    partition_sizes = Counter(sorted(partitions.values()))
    largest_partition = max(partition_sizes, key=partition_sizes.get)  # type: ignore

    for atom, partition in m.nodes(data=PARTITION):
        if partition_sizes[partition] == 1:
            # No need to artificially split singleton partitions.
            continue
        if partition != largest_partition:
            continue

        # Split largest partition.
        m_artificially_split = m.copy()
        nx.set_node_attributes(
            m_artificially_split,
            {atom: n_partitions},  # Partitions are zero-based.
            PARTITION,
        )
        m_artificially_refined = list(refine_partitions(m_artificially_split))[-1]

        yield m_artificially_refined


def get_discrete_partitionings(m: nx.Graph) -> Generator[nx.Graph, None, None]:
    """
    Build BFS refinement-tree and return its leaves (i.e., discrete partitionings).
    rtn = refinement-tree-node
    """
    rtn_queue = deque([m])

    while rtn_queue:
        rtn = rtn_queue.popleft()

        if partitioning_is_discrete(rtn):
            yield rtn
            continue

        rtn_queue.extend(list(get_refinement_tree_node_children(rtn)))


def get_canonical_molecule(ms: list[nx.Graph]) -> nx.Graph:
    m_canonical = None
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
    ms_discrete_partitionings = list(get_discrete_partitionings(m_refined))

    return get_canonical_molecule(ms_discrete_partitionings)
