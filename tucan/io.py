import networkx as nx
from pathlib import Path
from tucan.element_properties import ELEMENT_PROPS
from tucan.parser.parser import parse_tucan
from typing import Dict, List, Tuple


def graph_from_file(filepath: str) -> nx.Graph:
    """Instantiate a NetworkX graph from an MDL molfile.

    Parameters
    ----------
    filepath: str
        Path pointing to a file containing a v3000 MDL (now BIOVIA) connection
        table [1].

    Returns
    -------
    NetworkX Graph

    References
    ----------
    [1] https://en.wikipedia.org/wiki/Chemical_table_file
    """
    filepath = Path(filepath)
    if filepath.suffix != ".mol":
        raise IOError(f"The file must be in '.mol' format, not {filepath.suffix}.")
    filecontent = _read_file(filepath)

    return _graph_from_tokenized_lines(filecontent)


def graph_from_molfile_text(molfile: str) -> nx.Graph:
    lines = _split_into_tokenized_lines(molfile)
    return _graph_from_tokenized_lines(lines)


def graph_from_tucan(tucan: str) -> nx.Graph:
    return parse_tucan()


def _split_into_tokenized_lines(string: str) -> List[List[str]]:
    lines = [line.rstrip().split(" ") for line in string.splitlines()]
    return [[value for value in line if value != ""] for line in lines]


def _read_file(filepath: str) -> List[List[str]]:
    with open(filepath) as file:
        filecontent = file.read()
    return _split_into_tokenized_lines(filecontent)


def _graph_from_tokenized_lines(lines: List[List[str]]) -> nx.Graph:
    atom_props = _parse_atom_block_molfile3000(lines)
    bonds = _parse_bond_block_molfile3000(lines)

    graph = nx.Graph()
    graph.add_nodes_from(list(atom_props.keys()))
    nx.set_node_attributes(graph, atom_props)
    graph.add_edges_from(bonds)

    return graph


def _parse_atom_block_molfile3000(lines: List[List[str]]) -> Dict:
    atom_count = int(lines[5][3])
    atom_block_offset = 7

    atom_props = {
        int(line[2])
        - 1: _parse_atom_props(line)  # map zero-based label to property-dict
        for line in lines[atom_block_offset : atom_block_offset + atom_count]
    }
    assert len(atom_props) == atom_count, (
        f"Number of atoms {len(atom_props)} doesn't match atom-count specified in"
        f" header {atom_count}."
    )

    return atom_props


def _parse_atom_props(line: List[str]) -> Dict:
    element_symbol = line[3]
    isotope_mass = None
    if element_symbol == "D":
        element_symbol = "H"
        isotope_mass = 2
    elif element_symbol == "T":
        element_symbol = "H"
        isotope_mass = 3

    atom_props = {
        "element_symbol": element_symbol,
        "atomic_number": ELEMENT_PROPS[element_symbol]["atomic_number"],
        "partition": 0,
        "x_coord": float(line[4]),
        "y_coord": float(line[5]),
        "z_coord": float(line[6]),
    }

    optional_props = {
        "chg": [int(i.split("=")[1]) for i in line if "CHG" in i],
        "mass": [int(i.split("=")[1]) for i in line if "MASS" in i]
        if isotope_mass is None
        else [isotope_mass],
        "rad": [int(i.split("=")[1]) for i in line if "RAD" in i],
    }
    for key, val in optional_props.items():
        if val:
            atom_props[key] = val.pop()

    return atom_props


def _parse_bond_block_molfile3000(lines: List[List[str]]) -> List[Tuple[int, int]]:
    atom_count = int(lines[5][3])
    bond_count = int(lines[5][4])
    atom_block_offset = 7
    bond_block_offset = atom_block_offset + atom_count + 2

    bonds = [
        (int(line[4]) - 1, int(line[5]) - 1)
        for line in lines[bond_block_offset : bond_block_offset + bond_count]
    ]  # make bond-indices zero-based
    assert len(bonds) == bond_count, (
        f"Number of bonds {len(bonds)} doesn't match bond-count specified in header"
        f" {bond_count}."
    )

    return bonds
