from collections import Counter
from tucan.graph_utils import sort_molecule_by_attribute
import networkx as nx


def serialize_molecule(m):
    """Serialize a molecule."""
    m_sorted = sort_molecule_by_attribute(m, "atomic_number")
    serialization = _write_sum_formula(m_sorted)
    serialization += f"/{_write_edge_list(m_sorted)}"
    node_properties = _write_node_properties(m_sorted)
    serialization += f"/{node_properties}" if node_properties else ""

    return serialization


def _write_edge_list(m):
    sorted_edges = sorted([sorted(edge) for edge in m.edges()])
    edge_list_string = "".join(
        [f"({edge[0] + 1}-{edge[1] + 1})" for edge in sorted_edges]
    )

    return edge_list_string


def _write_node_properties(m):
    node_properties = ["chg", "mass", "rad"]
    node_property_string = ""
    for node in sorted(m.nodes(data=True)):
        label, props = node
        available_props = [
            f"{prop}={props[prop]}" for prop in node_properties if prop in props
        ]
        if not available_props:
            continue
        node_property_string += f"({label}:"
        node_property_string += f"{','.join(available_props)})"

    return node_property_string


def _write_sum_formula(m):
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
    sum_formula_string = ""
    carbon_count = element_counts.pop("C", None)
    if carbon_count:
        sum_formula_string += f"C{carbon_count}"
        hydrogen_count = element_counts.pop("H", None)
        if hydrogen_count:
            sum_formula_string += f"H{hydrogen_count}"
    for k, v in dict(sorted(element_counts.items())).items():
        sum_formula_string += f"{k}{v}"

    return sum_formula_string
