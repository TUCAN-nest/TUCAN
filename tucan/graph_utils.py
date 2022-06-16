import networkx as nx
import random
from tucan.io import graph_from_moldata


def sort_molecule_by_attribute(m, attribute):
    """Sort atoms by attribute."""
    attr_sequence = [attribute_sequence(atom, m, attribute) for atom in m]
    attr_with_labels = [
        (i, j) for i, j in zip(attr_sequence, list(m.nodes))
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

    element_symbols_permu = [e for (i, e) in m_permu.nodes(data="element_symbol")]
    return graph_from_moldata(
        list(m_permu.nodes), element_symbols_permu, list(m_permu.edges)
    )


def _permute_molecule(m):
    labels = list(m.nodes)
    permuted_labels = list(labels)  # shallow copy
    random.shuffle(permuted_labels)
    return relabel_molecule(m, permuted_labels, labels)
