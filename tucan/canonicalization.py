from tucan.graph_attributes import PARTITION
from tucan.graph_utils import attribute_sequence
import networkx as nx
from igraph import Graph as iGraph
from typing import Iterator


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
    return max(nx.get_node_attributes(m, PARTITION).values())


def refine_partitions(m: nx.Graph) -> Iterator[nx.Graph]:
    n_current_partitions = get_number_of_partitions(m)

    if n_current_partitions == m.number_of_nodes() - 1:
        # partitions are discrete (i.e., each node in a separate partition)
        return m

    m_refined = partition_molecule_by_attribute(m, PARTITION)
    if get_number_of_partitions(m_refined) == n_current_partitions:
        # no refinement possible
        return m

    yield m_refined

    yield from refine_partitions(m_refined)


def assign_canonical_labels(m: nx.Graph) -> dict[int, int]:
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
    partitions = m_igraph.vs[PARTITION]
    canonical_labels = m_igraph.canonical_permutation(color=partitions)

    return dict(zip(old_labels, canonical_labels))


def canonicalize_molecule(m: nx.Graph) -> nx.Graph:
    m_partitioned_by_invariant_code = partition_molecule_by_attribute(
        m, "invariant_code"
    )
    m_refined = list(refine_partitions(m_partitioned_by_invariant_code))
    m_partitioned = m_refined[-1] if m_refined else m_partitioned_by_invariant_code

    canonical_labels = assign_canonical_labels(m_partitioned)

    return nx.relabel_nodes(m_partitioned, canonical_labels, copy=True)
