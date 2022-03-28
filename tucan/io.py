import networkx as nx
import random
from tucan.element_properties import ELEMENT_PROPS
from tucan.canonicalization import (
    _cycle_memberships,
    _add_invariant_code,
    _relabel_molecule,
)
from rdkit import Chem
from typing import List, Tuple


def graph_from_file(filepath):
    with open(filepath) as f:
        filecontent = f.read()
    if filepath.suffix == ".mol":
        element_symbols, bonds = _parse_molfile3000(filecontent)
    elif filepath.suffix == ".col":
        element_symbols, bonds = _parse_dimacs(filecontent)
    else:
        raise IOError("Invalid file format, must be one of {.mol, .col}.")
    return _graph_from_moldata(element_symbols, bonds)


def graph_from_smiles(smiles: str):
    molfile = _molfile3000_from_smiles(smiles)
    element_symbols, bonds = _parse_molfile3000(molfile)
    return _graph_from_moldata(element_symbols, bonds)


def _molfile3000_from_smiles(smiles: str):
    m = Chem.MolFromSmiles(smiles, sanitize=False)
    return Chem.MolToMolBlock(m, forceV3000=True, includeStereo=False, kekulize=False)


def _graph_from_moldata(element_symbols: List[str], bonds: List[Tuple[int]]):
    """Instantiate a NetworkX graph from molecular data.

    Parameters
    ----------
    element_symbols: List[str]
        Element symbols associated with the atoms.
    bonds: List[Tuple[int, int]]
        Bonds between atoms.
    """
    atomic_numbers = [ELEMENT_PROPS[s]["atomic_number"] for s in element_symbols]
    node_labels = range(len(element_symbols))
    graph = nx.Graph()
    graph.add_nodes_from(node_labels)
    graph.add_edges_from(bonds)
    nx.set_node_attributes(
        graph, dict(zip(node_labels, element_symbols)), "element_symbol"
    )
    nx.set_node_attributes(
        graph, dict(zip(node_labels, atomic_numbers)), "atomic_number"
    )
    nx.set_node_attributes(graph, _cycle_memberships(graph), "cycle_membership")
    nx.set_node_attributes(graph, 0, "partition")
    _add_invariant_code(graph)
    return graph


def _parse_molfile3000(filecontent: str):
    lines = [l.rstrip().split(" ") for l in filecontent.splitlines()]
    lines = [[value for value in line if value != ""] for line in lines]
    atom_count = int(lines[5][3])
    bond_count = int(lines[5][4])
    atom_block_offset = 7
    bond_block_offset = atom_block_offset + atom_count + 2
    element_symbols = [
        l[3] for l in lines[atom_block_offset : atom_block_offset + atom_count]
    ]
    assert (
        len(element_symbols) == atom_count
    ), f"Number of atoms {len(element_symbols)} doesn't match atom-count specified in header {atom_count}."
    bonds = [
        (int(l[4]) - 1, int(l[5]) - 1)
        for l in lines[bond_block_offset : bond_block_offset + bond_count]
    ]  # make bond-indices zero-based
    assert (
        len(bonds) == bond_count
    ), f"Number of bonds {len(bonds)} doesn't match bond-count specified in header {bond_count}."
    return element_symbols, bonds


def _parse_molfile2000(filecontent: str):
    lines = [
        [l[i : i + 3].strip(" ") for i in range(0, len(l), 3)]
        for l in filecontent.splitlines()
    ]
    atom_count = int(lines[3][0])
    bond_count = int(lines[3][1])
    atom_block_offset = 4
    bond_block_offset = atom_block_offset + atom_count
    element_symbols = [
        l[10] for l in lines[atom_block_offset : atom_block_offset + atom_count]
    ]
    assert (
        len(element_symbols) == atom_count
    ), f"Number of atoms {len(element_symbols)} doesn't match atom-count specified in header {atom_count}."
    bonds = [
        (int(l[0]) - 1, int(l[1]) - 1)
        for l in lines[bond_block_offset : bond_block_offset + bond_count]
    ]  # make bond-indices zero-based
    assert (
        len(bonds) == bond_count
    ), f"Number of bonds {len(bonds)} doesn't match bond-count specified in header {bond_count}."
    return element_symbols, bonds


def _parse_dimacs(filecontent: str):
    lines = [l.rstrip().split(" ") for l in filecontent.splitlines()]
    lines = [l for l in lines if l[0] in ["p", "e"]]
    atom_count = int(lines[0][2])
    element_symbols = ["C"] * atom_count
    bonds = [
        (int(l[1]) - 1, int(l[2]) - 1) for l in lines[1:]
    ]  # make bond-indices zero-based
    return element_symbols, bonds


def permute_molecule(m, random_seed=1.0):
    """Randomly permute the atom-labels of a molecular graph.

    Parameters
    ----------
    random_seed: float
        In [0.0, 1.0).
    """
    labels = m.nodes()
    permuted_labels = list(range(m.number_of_nodes()))
    # Enforce permutation for graphs with at least 2 edges that aren't fully connected (i.e., complete).
    enforce_permutation = m.number_of_edges() > 1 and nx.density(m) != 1
    random.seed(
        random_seed
    )  # subsequent calls of random.shuffle(x[, random]) will now use fixed sequence of values for `random` parameter
    random.shuffle(permuted_labels)
    m_permu = _relabel_molecule(m, permuted_labels, labels)
    if enforce_permutation:
        while m.edges == m_permu.edges:
            random.shuffle(permuted_labels)
            m_permu = _relabel_molecule(m, permuted_labels, labels)
    return m_permu
