from __future__ import annotations
import networkx as nx
import random
from typing import Any, NamedTuple

from tucan.graph_attributes import (
    ATOMIC_NUMBER,
    INVARIANT_CODE,
    MASS,
    RAD,
)


def graph_from_molecule(
    atom_attrs: dict[int, dict[str, Any]],
    bond_attrs: dict[tuple[int, int], dict[str, int]],
) -> nx.Graph:
    invariant_code_definitions = [
        InvariantCodeDefinition(ATOMIC_NUMBER),
        InvariantCodeDefinition(MASS, 0),
        InvariantCodeDefinition(RAD, 0),
    ]
    _add_invariant_code(atom_attrs, invariant_code_definitions)

    graph = nx.Graph()
    graph.add_nodes_from(list(atom_attrs.keys()))
    nx.set_node_attributes(graph, atom_attrs)
    graph.add_edges_from(list(bond_attrs.keys()))
    nx.set_edge_attributes(graph, bond_attrs)

    return nx.convert_node_labels_to_integers(graph)


class InvariantCodeDefinition(NamedTuple):
    key: str
    default_value: Any = None


def _add_invariant_code(
    atom_attrs: dict[int, dict[str, Any]],
    invariant_code_definitions: list[InvariantCodeDefinition],
) -> None:
    for atom, attrs in atom_attrs.items():
        invariant_code = tuple(
            (
                attrs[icd.key]
                if (default_value := icd.default_value) is None
                else attrs.get(icd.key, default_value)
            )
            for icd in invariant_code_definitions
        )
        atom_attrs[atom].update({INVARIANT_CODE: invariant_code})


def get_attribute_sequences(
    attributes: dict[int, Any], neighbors: dict[int, tuple[int]]
) -> tuple[tuple[Any, ...], ...]:
    return tuple(
        (attribute, *neighbor_attrs)
        for node, attribute in attributes.items()
        for neighbor_attrs in [[attributes[neighbor] for neighbor in neighbors[node]]]
        for _ in [neighbor_attrs.sort()]
    )


def sort_molecule_by_attribute(m: nx.Graph, attribute: str) -> nx.Graph:
    """Sort atoms by attribute."""
    attributes = dict(nx.get_node_attributes(m, attribute))
    neighbors = {node: tuple(m[node]) for node in m}

    attr_with_labels = zip(
        get_attribute_sequences(attributes, neighbors), attributes.keys(), strict=True
    )  # [(A, 0), (C, 1), (B, 2)]
    sorted_attr, labels_sorted_by_attr = zip(
        *sorted(attr_with_labels), strict=True
    )  # (A, B, C), (0, 2, 1)

    return nx.relabel_nodes(
        m,
        dict(zip(labels_sorted_by_attr, list(range(m.number_of_nodes())), strict=True)),
        copy=True,
    )


def permute_molecule(m: nx.Graph, random_seed: float = 1.0) -> nx.Graph:
    """Randomly permute the atom-labels of a molecular graph.

    Parameters
    ----------
    random_seed: float
        In [0.0, 1.0).
    """
    random.seed(
        random_seed
    )  # subsequent calls of random.shuffle(x[, random]) will now use fixed sequence of values for `random` parameter

    m_permu = _permute_molecule(m)

    # Enforce permutation for graphs with at least 2 edges that aren't fully connected (i.e., complete).
    enforce_permutation = m.number_of_edges() > 1 and nx.density(m) != 1
    if enforce_permutation:
        while m.edges == m_permu.edges:
            m_permu = _permute_molecule(m)

    return m_permu


def _permute_molecule(m: nx.Graph) -> nx.Graph:
    labels = list(m.nodes)
    permuted_labels = list(labels)  # shallow copy
    random.shuffle(permuted_labels)
    m_relabeled = nx.relabel_nodes(
        m, dict(zip(permuted_labels, labels, strict=True)), copy=True
    )

    return _sort_molecule_by_label(m_relabeled)


def _sort_molecule_by_label(m: nx.Graph) -> nx.Graph:
    """Sort molecule by label.

    Ensure that the graph's node iteration order is identical to the label order.
    In a NetworkX graph, the iteration order of the nodes depends on the initial
    insertion order of the nodes. There's a crucial difference to `nx.relabel_nodes()`.
    The latter only changes the labels, without changing the iteration order.

    Original nodes:
        label | attribute
        -----------------
        0     | A
        1     | B
        2     | C

    Nodes relabeled with `nx.relabel_nodes()`:
        label | attribute
        -----------------
        1     | A
        2     | B
        0     | C

    In contrast, the present function changes the labels _and_ iteration order:
        label | attribute
        -----------------
        0     | C
        1     | A
        2     | B
    """
    nodes_sorted_by_label = sorted(list(m.nodes(data=True)))

    m_sorted_by_label = nx.Graph()
    m_sorted_by_label.add_nodes_from(nodes_sorted_by_label)
    m_sorted_by_label.add_edges_from(m.edges(data=True))

    return m_sorted_by_label
