from collections import Counter
from tucan.graph_utils import sort_molecule_by_attribute
from operator import gt, lt, eq
import networkx as nx


def serialize_molecule(m):
    """Serialize a molecule."""
    m_sorted = sort_molecule_by_attribute(_assign_final_labels(m), "atomic_number")
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
    m,
    traversal_priorities=[lt, gt, eq],
):
    """Re-label nodes such that each node is connected to neighbors with the
    smallest possible labels.
    This is not part of (and not required for) the canonicalization.
    The re-labeling is for cosmetic purposes."""
    partitions = m.nodes.data("partition")
    labels_by_partition = _labels_by_partition(m)
    atom_queue = [0]
    final_labels = {}
    nx.set_node_attributes(m, False, "explored")

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

        for n in neighbor_traversal_order:
            atom_queue.insert(0, n)

    nx.set_node_attributes(m, False, "explored")

    return nx.relabel_nodes(m, final_labels, copy=True)


def _labels_by_partition(m):
    """Create dictionary of partitions to node labels."""
    partitions = set(sorted([v for _, v in m.nodes.data("partition")]))
    labels_by_partition = {p: set() for p in partitions}
    for a in m:
        labels_by_partition[m.nodes[a]["partition"]].add(a)
    labels_by_partition.update(
        (k, sorted(list(v), reverse=True)) for k, v in labels_by_partition.items()
    )
    return labels_by_partition
