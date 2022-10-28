import networkx as nx
import random


def sort_molecule_by_attribute(m, attribute):
    """Sort atoms by attribute."""
    attr_with_labels = [
        (attribute_sequence(atom, m, attribute), atom) for atom in m
    ]  # [(A, 0), (C, 1), (B, 2)]
    sorted_attr, labels_sorted_by_attr = zip(
        *sorted(attr_with_labels)
    )  # (A, B, C), (0, 2, 1)
    return relabel_molecule(m, labels_sorted_by_attr, list(range(m.number_of_nodes())))


def attribute_sequence(atom, m, attribute):
    attr_atom = m.nodes[atom][attribute]
    attr_neighbors = sorted(
        [m.nodes[n][attribute] for n in m.neighbors(atom)], reverse=True
    )
    return [attr_atom] + attr_neighbors


def relabel_molecule(m, old_labels, new_labels):
    """Relabel the atoms of a molecular graph."""
    return nx.relabel_nodes(m, dict(zip(old_labels, new_labels)), copy=True)


def permute_molecule(m, random_seed=1.0):
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


def _permute_molecule(m):
    labels = list(m.nodes)
    permuted_labels = list(labels)  # shallow copy
    random.shuffle(permuted_labels)
    m_relabeled = relabel_molecule(m, permuted_labels, labels)
    return _sort_molecule_by_label(m_relabeled)


def _sort_molecule_by_label(m):
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
