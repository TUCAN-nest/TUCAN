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
