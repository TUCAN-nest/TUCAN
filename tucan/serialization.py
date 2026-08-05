from collections import Counter

from tucan.graph_attributes import (
    ATOMIC_NUMBER,
    ELEMENT_SYMBOL,
    MASS,
    RAD,
)
from tucan.graph_utils import sort_molecule_by_attribute
from typing import Final
import networkx as nx


def serialize_molecule(m: nx.Graph) -> str:
    """Serialize a molecule."""
    m_sorted = sort_molecule_by_attribute(m, ATOMIC_NUMBER)
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


_SERIALIZER_NODE_ATTRIBUTE_MAPPING: Final[dict[str, str]] = {
    MASS: "mass",
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
