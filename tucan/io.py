import networkx as nx
from tucan.element_properties import ELEMENT_PROPS
from typing import List, Tuple


def graph_from_file(filepath):
    with open(filepath) as f:
        filecontent = f.read()
    if filepath.suffix == ".mol":
        node_labels, element_symbols, bonds = _parse_molfile3000(filecontent)
    elif filepath.suffix == ".col":
        node_labels, element_symbols, bonds = _parse_dimacs(filecontent)
    else:
        raise IOError("Invalid file format, must be one of {.mol, .col}.")
    return graph_from_moldata(node_labels, element_symbols, bonds)


def graph_from_moldata(
    node_labels: List[int], element_symbols: List[str], bonds: List[Tuple[int]]
):
    """Instantiate a NetworkX graph from molecular data.

    Parameters
    ----------
    node_labels: List[int]
        Node labels associated with the atoms.
    element_symbols: List[str]
        Element symbols associated with the atoms.
    bonds: List[Tuple[int, int]]
        Bonds between atoms.
    """
    atomic_numbers = [ELEMENT_PROPS[s]["atomic_number"] for s in element_symbols]
    graph = nx.Graph()
    graph.add_nodes_from(node_labels)
    graph.add_edges_from(bonds)
    nx.set_node_attributes(
        graph, dict(zip(node_labels, element_symbols)), "element_symbol"
    )
    nx.set_node_attributes(
        graph, dict(zip(node_labels, atomic_numbers)), "atomic_number"
    )
    nx.set_node_attributes(graph, 0, "partition")
    return graph


def _parse_molfile3000(filecontent: str):
    lines = [l.rstrip().split(" ") for l in filecontent.splitlines()]
    lines = [[value for value in line if value != ""] for line in lines]
    atom_count = int(lines[5][3])
    bond_count = int(lines[5][4])
    atom_block_offset = 7
    bond_block_offset = atom_block_offset + atom_count + 2
    node_labels = [
        int(l[2]) - 1 for l in lines[atom_block_offset : atom_block_offset + atom_count]
    ]  # convert from 1-based to 0-based
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
    return node_labels, element_symbols, bonds


def _parse_dimacs(filecontent: str):
    lines = [l.rstrip().split(" ") for l in filecontent.splitlines()]
    lines = [l for l in lines if l[0] in ["p", "e"]]
    atom_count = int(lines[0][2])
    node_labels = range(atom_count)
    element_symbols = ["C"] * atom_count
    bonds = [
        (int(l[1]) - 1, int(l[2]) - 1) for l in lines[1:]
    ]  # make bond-indices zero-based
    return node_labels, element_symbols, bonds
