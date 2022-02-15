from rdkit.Chem import rdmolfiles
from collections import Counter
import networkx as nx
import random


def relabel_molecule(m, old_labels, new_labels):
    """Relabel the atoms of a molecular graph."""
    m_relabeled = nx.relabel_nodes(m, dict(zip(old_labels, new_labels)))
    # In the NetworkX Graph datastructure, the relabeled nodes don't occur in
    # increasing order yet. This is why we change the node order now.
    m_sorted = nx.Graph()
    m_sorted.add_nodes_from(sorted(m_relabeled.nodes(data=True)))
    m_sorted.add_edges_from(m_relabeled.edges(data=True))
    return m_sorted


def permute_molecule(m, random_seed=1.0):
    """Randomly permute the atom-labels of a molecular graph.

    Parameters
    ----------
    random_seed: float
        In [0.0, 1.0).
    """
    idcs = m.nodes()
    permuted_idcs = list(range(m.number_of_nodes()))
    # Enforce permutation for graphs with at least 2 edges that aren't fully connected (i.e., complete).
    enforce_permutation = m.number_of_edges() > 1 and nx.density(m) != 1
    random.seed(
        random_seed
    )  # subsequent calls of random.shuffle(x[, random]) will now use fixed sequence of values for `random` parameter
    random.shuffle(permuted_idcs)
    m_permu = relabel_molecule(m, permuted_idcs, idcs)
    if enforce_permutation:
        while m.edges == m_permu.edges:
            random.shuffle(permuted_idcs)
            m_permu = relabel_molecule(m, permuted_idcs, idcs)
    return m_permu


def serialize_molecule(m):
    """Serialize a molecule."""
    serialization = write_sum_formula(m)
    for edge in sorted([sorted(edge) for edge in m.edges()]):
        serialization += f"/{edge[0]}-{edge[1]}"
    return serialization


def write_sum_formula(m):
    """Write sum formula of a molecule.

    Elements occur in the following order: C, H, other elements in alphabetic order.
    """
    element_counts = Counter(nx.get_node_attributes(m, "element_symbol").values())
    element_counts = {
        k: (v if v > 1 else "") for k, v in element_counts.items()
    }  # remove counts of 1 since those are implicit in sum formula
    sum_formula = ""
    for element in ["C", "H"]:
        count = element_counts.pop(element, None)
        if count:
            sum_formula += f"{element}{count}"
    for k, v in dict(sorted(element_counts.items())).items():
        sum_formula += f"{k}{v}"
    return sum_formula


def mdl_2000_to_3000(molfile_path, removeHs=False):
    v2000 = rdmolfiles.MolFromMolFile(molfile_path, removeHs=removeHs)
    rdmolfiles.MolToMolFile(v2000, molfile_path, forceV3000=True)
