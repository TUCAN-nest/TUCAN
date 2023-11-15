from collections import Counter, deque

from tucan.graph_attributes import (
    ATOMIC_NUMBER,
    ELEMENT_SYMBOL,
    PARTITION,
    RAD,
)
from tucan.graph_utils import sort_molecule_by_attribute
from operator import gt, lt, eq
from typing import Callable
import networkx as nx


def serialize_molecule(m: nx.Graph) -> str:
    """Serialize a molecule."""
    m_sorted = sort_molecule_by_attribute(_assign_final_labels(m), ATOMIC_NUMBER)
    serialization = _write_sum_formula(m_sorted)
    serialization += f"/{_write_edge_list(m_sorted)}"
    node_attributes = _write_node_attributes(m_sorted)
    serialization += f"/{node_attributes}" if node_attributes else ""

    return serialization


def _write_edge_list(m: nx.Graph) -> str:
    sorted_edges = sorted([sorted(edge) for edge in m.edges()])
    edge_list_string = "".join(
        [f"({edge[0] + 1}-{edge[1] + 1})" for edge in sorted_edges]
    )

    return edge_list_string


_SERIALIZER_NODE_ATTRIBUTE_MAPPING = {
    "mass": "mass",
    RAD: "rad",
}


def _write_node_attributes(m: nx.Graph) -> str:
    node_attribute_string = ""
    for label, attrs in sorted(m.nodes(data=True)):
        available_attrs = [
            f"{_SERIALIZER_NODE_ATTRIBUTE_MAPPING[attr]}={attrs[attr]}"
            for attr in _SERIALIZER_NODE_ATTRIBUTE_MAPPING
            if attr in attrs
        ]
        if not available_attrs:
            continue
        node_attribute_string += f"({label + 1}:"
        node_attribute_string += f"{','.join(available_attrs)})"

    return node_attribute_string


def _write_sum_formula(m: nx.Graph) -> str:
    """Write the sum formula of a molecule.

    In the sum formula the elements are ordered according to Hill system [1]:
    1. C
    2. H
    3. symbols of remaining elements in alphabetic order (including H in other
    than carbon compounds)

    References
    ----------
    [1] doi:10.1021/ja02046a005
    """
    element_counts = Counter(nx.get_node_attributes(m, ELEMENT_SYMBOL).values())
    sum_formula_string = ""
    carbon_count = element_counts.pop("C", None)
    if carbon_count:
        sum_formula_string += f"C{carbon_count}" if carbon_count > 1 else "C"
        hydrogen_count = element_counts.pop("H", None)
        if hydrogen_count:
            sum_formula_string += f"H{hydrogen_count}" if hydrogen_count > 1 else "H"
    for k, v in dict(sorted(element_counts.items())).items():
        sum_formula_string += f"{k}{v}" if v > 1 else k

    return sum_formula_string


def _assign_final_labels(
    m: nx.Graph,
    traversal_priorities: tuple[Callable, Callable, Callable] = (lt, gt, eq),
) -> nx.Graph:
    """Re-label nodes such that each node is connected to neighbors with the
    smallest possible labels.
    This is not part of (and not required for) the canonicalization.
    The re-labeling is for cosmetic purposes."""
    partitions = m.nodes.data(PARTITION)
    labels_by_partition = _labels_by_partition(m)
    final_labels = {}
    nx.set_node_attributes(m, False, "explored")

    # outer loop iterates over all fragments of the graph (= graph components),
    # starting with the lowest unexplored node label
    while unexplored := sorted([k for k, v in m.nodes(data="explored") if not v]):
        atom_queue = deque([unexplored[0]])

        # inner loop reaches out to all atoms in a fragment
        while atom_queue:
            a = atom_queue.pop()
            if m.nodes[a]["explored"]:
                continue
            a_final = labels_by_partition[
                partitions[a]
            ].pop()  # assign smallest label available in this partition
            final_labels[a] = a_final

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

            atom_queue.extendleft(neighbor_traversal_order)

    assert len(final_labels) == len(m.nodes)

    nx.set_node_attributes(m, False, "explored")
    return nx.relabel_nodes(m, final_labels, copy=True)


def _labels_by_partition(m: nx.Graph) -> dict[int, list[int]]:
    """Create dictionary of partitions to node labels."""
    partitions = set(sorted([v for _, v in m.nodes.data(PARTITION)]))
    labels_by_partition: dict[int, list[int]] = {p: [] for p in partitions}
    for a in m:
        labels_by_partition[m.nodes[a][PARTITION]].append(a)
    labels_by_partition.update(
        (k, sorted(list(v), reverse=True)) for k, v in labels_by_partition.items()
    )
    return labels_by_partition
