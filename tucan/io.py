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
    """Instantiate a NetworkX graph from a TUCAN string.

    Parameters
    ----------
    tucan: str
        TUCAN string to be deserialized.

    Returns
    -------
    NetworkX Graph
    """
    return parse_tucan()


def _split_into_tokenized_lines(string: str) -> List[List[str]]:
    lines = string.splitlines()
    lines = _concat_lines_with_dash(lines)
    splitted_lines = [line.rstrip().split(" ") for line in lines]
    return [[value for value in line if value != ""] for line in splitted_lines]


def _concat_lines_with_dash(lines: List[str]) -> List[str]:
    result = list()

    while (length := len(lines)) > 1:
        curr_line = lines.pop()
        prev_line = lines[length - 2]

        if prev_line.endswith("-") and prev_line.startswith("M  V30 "):
            if curr_line.startswith("M  V30 "):
                lines[length - 2] = prev_line[0:-1] + curr_line[7:]
            else:
                raise MolfileParserException(
                    f'Invalid concatenation of lines "{prev_line}" and "{curr_line}"'
                )
        else:
            result.append(curr_line)

    result.append(lines[0])
    result.reverse()
    return result


def _read_file(filepath: str) -> List[List[str]]:
    with open(filepath) as file:
        filecontent = file.read()
    return _split_into_tokenized_lines(filecontent)


def _graph_from_tokenized_lines(lines: List[List[str]]) -> nx.Graph:
    _validate_molfile_version(lines, "V3000")
    _validate_counts_line(lines)

    atom_props = _parse_atom_block_molfile3000(lines)
    bonds = _parse_bond_block_molfile3000(lines)

    _validate_bond_indices(bonds, atom_props)

    graph = nx.Graph()
    graph.add_nodes_from(list(atom_props.keys()))
    nx.set_node_attributes(graph, atom_props)
    graph.add_edges_from(bonds)

    return graph


def _validate_molfile_version(lines: List[List[str]], expected_version: str):
    if (version := lines[3][6]) != expected_version:
        raise MolfileParserException(
            f'Invalid Molfile version: Expected "{expected_version}", found "{version}"'
        )


def _validate_counts_line(lines: List[List[str]]):
    if lines[5][2] != "COUNTS" or len(lines[5]) < 5:
        badline = " ".join(lines[5])
        raise MolfileParserException(f'Bad counts line: "{badline}"')


def _parse_atom_block_molfile3000(lines: List[List[str]]) -> Dict:
    atom_count = int(lines[5][3])
    atom_block_offset = 7

    if (begin_atom_str := " ".join(lines[atom_block_offset - 1][2:])) != "BEGIN ATOM":
        raise MolfileParserException(
            f'Expected "BEGIN ATOM" in line {atom_block_offset}, found "{begin_atom_str}"'
        )
    if (
        end_atom_str := " ".join(lines[atom_block_offset + atom_count][2:])
    ) != "END ATOM":
        raise MolfileParserException(
            f'Expected "END ATOM" in line {atom_block_offset + atom_count + 1}, found "{end_atom_str}"'
        )

    atom_props = {
        int(line[2])
        - 1: _parse_atom_props(line)  # map zero-based label to property-dict
        for line in lines[atom_block_offset : atom_block_offset + atom_count]
    }

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

    # bond block is optional
    if bond_count == 0:
        return []

    atom_block_offset = 7
    bond_block_offset = atom_block_offset + atom_count + 2

    if (begin_bond_str := " ".join(lines[bond_block_offset - 1][2:])) != "BEGIN BOND":
        raise MolfileParserException(
            f'Expected "BEGIN BOND" in line {bond_block_offset}, found "{begin_bond_str}"'
        )
    if (
        end_bond_str := " ".join(lines[bond_block_offset + bond_count][2:])
    ) != "END BOND":
        raise MolfileParserException(
            f'Expected "END BOND" in line {bond_block_offset + bond_count + 1}, found "{end_bond_str}"'
        )

    bonds = [
        (int(line[4]) - 1, int(line[5]) - 1)
        for line in lines[bond_block_offset : bond_block_offset + bond_count]
    ]  # make bond-indices zero-based

    return bonds


def _validate_bond_indices(bonds: List[Tuple[int, int]], atom_props: Dict):
    for bond in bonds:
        _validate_atom_index(bond[0], atom_props)
        _validate_atom_index(bond[1], atom_props)


def _validate_atom_index(index: int, atom_props: Dict):
    if index not in atom_props:
        raise MolfileParserException(f"Unknown atom index {index + 1} in bond")


class MolfileParserException(Exception):
    pass
