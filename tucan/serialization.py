from collections import Counter
from tucan.graph_utils import sort_molecule_by_attribute
import networkx as nx


def serialize_molecule(m):
    """Serialize a molecule."""
    serialization = _sum_formula(m)
    m_sorted_by_atomic_number = sort_molecule_by_attribute(m, "atomic_number")
    for edge in _edge_list(m_sorted_by_atomic_number):
        serialization += f"/{edge[0] + 1}-{edge[1] + 1}"
    return serialization


def _edge_list(m):
    return sorted([sorted(edge) for edge in m.edges()])


def _sum_formula(m):
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
    element_counts = Counter(nx.get_node_attributes(m, "element_symbol").values())
    element_counts = {
        k: (v if v > 1 else "") for k, v in element_counts.items()
    }  # remove counts of 1 since those are implicit in sum formula
    sum_formula = ""
    carbon_count = element_counts.pop("C", None)
    if carbon_count:
        sum_formula += f"C{carbon_count}"
        hydrogen_count = element_counts.pop("H", None)
        if hydrogen_count:
            sum_formula += f"H{hydrogen_count}"
    for k, v in dict(sorted(element_counts.items())).items():
        sum_formula += f"{k}{v}"
    return sum_formula
