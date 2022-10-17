import networkx as nx
from collections import deque
from pathlib import Path
from tucan.element_properties import ELEMENT_PROPS


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


def _split_into_tokenized_lines(string: str) -> list[list[str]]:
    lines = string.splitlines()
    lines = _concat_lines_with_dash(lines)
    split_lines = [line.rstrip().split(" ") for line in lines]
    return [[value for value in line if value != ""] for line in split_lines]


def _concat_lines_with_dash(lines: list[str]) -> list[str]:
    final_lines = []
    lines = deque(lines)

    while lines:
        curr_line = lines.popleft()

        if not lines:
            final_lines.append(curr_line)
            break

        if not (curr_line.startswith("M  V30 ") and curr_line.endswith("-")):
            final_lines.append(curr_line)
            continue

        next_line = lines.popleft()
        if not next_line.startswith("M  V30 "):
            raise MolfileParserException(
                f'Invalid concatenation of lines "{curr_line}" and "{next_line}"'
            )

        lines.appendleft(curr_line[0:-1] + next_line[7:])

    return final_lines


def _read_file(filepath: str) -> list[list[str]]:
    with open(filepath) as file:
        filecontent = file.read()
    return _split_into_tokenized_lines(filecontent)


def _graph_from_tokenized_lines(lines: list[list[str]]) -> nx.Graph:
    _validate_molfile_version(lines, "V3000")
    _validate_counts_line(lines)

    atom_props, star_atoms = _parse_atom_block_molfile3000(lines)
    bond_props = _parse_bond_block_molfile3000(lines, star_atoms)

    _validate_bond_indices(bond_props, atom_props)

    graph = nx.Graph()
    graph.add_nodes_from(list(atom_props.keys()))
    nx.set_node_attributes(graph, atom_props)
    graph.add_edges_from(list(bond_props.keys()))
    nx.set_edge_attributes(graph, bond_props)

    mapping = {old_index: new_index for new_index, old_index in enumerate(graph.nodes)}
    nx.relabel_nodes(graph, mapping, copy=False)

    return graph


def _validate_molfile_version(lines: list[list[str]], expected_version: str):
    if (version := lines[3][-1]) != expected_version:
        raise MolfileParserException(
            f'Invalid Molfile version: Expected "{expected_version}", found "{version}"'
        )


def _validate_counts_line(lines: list[list[str]]):
    if lines[5][2] != "COUNTS" or len(lines[5]) < 5:
        badline = " ".join(lines[5])
        raise MolfileParserException(f'Bad counts line: "{badline}"')


def _parse_atom_block_molfile3000(lines: list[list[str]]) -> tuple[dict, list[int]]:
    atom_count = int(lines[5][3])
    atom_block_offset = 7

    if (begin_atom_str := " ".join(lines[atom_block_offset - 1][2:])) != "BEGIN ATOM":
        raise MolfileParserException(
            f'Expected "BEGIN ATOM" on line {atom_block_offset}, found "{begin_atom_str}"'
        )
    if (
        end_atom_str := " ".join(lines[atom_block_offset + atom_count][2:])
    ) != "END ATOM":
        raise MolfileParserException(
            f'Expected "END ATOM" on line {atom_block_offset + atom_count + 1}, found "{end_atom_str}"'
        )

    atom_props = {}
    star_atoms = []
    for line in lines[atom_block_offset : atom_block_offset + atom_count]:
        atom_index = int(line[2]) - 1
        props, is_star_atom = _parse_atom_props(line)
        if is_star_atom:
            star_atoms.append(atom_index)
            continue
        atom_props[atom_index] = props

    return atom_props, star_atoms


def _parse_atom_props(line: list[str]) -> tuple[dict, bool]:
    element_symbol = line[3]
    if element_symbol == "*":
        return None, True

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

    return atom_props, False


def _parse_bond_block_molfile3000(
    lines: list[list[str]], star_atoms: list[int]
) -> dict[tuple[int, int]]:
    atom_count = int(lines[5][3])
    bond_count = int(lines[5][4])

    # bond block is optional
    if bond_count == 0:
        return {}

    atom_block_offset = 7
    bond_block_offset = atom_block_offset + atom_count + 2

    if (begin_bond_str := " ".join(lines[bond_block_offset - 1][2:])) != "BEGIN BOND":
        raise MolfileParserException(
            f'Expected "BEGIN BOND" on line {bond_block_offset}, found "{begin_bond_str}"'
        )
    if (
        end_bond_str := " ".join(lines[bond_block_offset + bond_count][2:])
    ) != "END BOND":
        raise MolfileParserException(
            f'Expected "END BOND" on line {bond_block_offset + bond_count + 1}, found "{end_bond_str}"'
        )

    bonds = {}
    for line in lines[bond_block_offset : bond_block_offset + bond_count]:
        atom1_index = int(line[4]) - 1
        atom2_index = int(line[5]) - 1
        bond_props = _parse_bond_props(line)

        atom1_is_star = atom1_index in star_atoms
        atom2_is_star = atom2_index in star_atoms
        if atom1_is_star and atom2_is_star:
            raise MolfileParserException(
                f'Two "star" atoms (index {atom1_index + 1} and {atom2_index + 1}) may not be connected'
            )
        elif atom1_is_star:
            bond_tuples = _parse_bond_line_with_star_atom(line, atom2_index)
        elif atom2_is_star:
            bond_tuples = _parse_bond_line_with_star_atom(line, atom1_index)
        else:
            bond_tuples = [(atom1_index, atom2_index)]

        for t in bond_tuples:
            bonds[t] = bond_props.copy()

    return bonds


def _parse_bond_props(line: list[str]) -> dict:
    return {"bond_type": int(line[3])}


def _parse_bond_line_with_star_atom(
    line: list[str], start_atom_index: int
) -> list[tuple[int, int]]:
    for index, token in enumerate(line):
        if token.startswith("ENDPTS=("):
            endpts_index = index
            number_of_endpts = int(token.split("=(")[1])
            break
    # silently ignore everything that has no ENDPTS (e.g. use of star atoms in polymers)
    else:
        return []

    endpoints = [
        int(token) for token in line[endpts_index + 1 : endpts_index + number_of_endpts]
    ]
    if not line[endpts_index + number_of_endpts].endswith(")"):
        raise MolfileParserException("Expected end of ENDPTS block")
    endpoints.append(int(line[endpts_index + number_of_endpts][:-1]))

    return [(start_atom_index, end_atom_index - 1) for end_atom_index in endpoints]


def _validate_bond_indices(bonds: dict[tuple[int, int]], atom_props: dict):
    for bond in bonds.keys():
        _validate_atom_index(bond[0], atom_props)
        _validate_atom_index(bond[1], atom_props)


def _validate_atom_index(index: int, atom_props: dict):
    if index not in atom_props:
        raise MolfileParserException(f"Unknown atom index {index + 1} in bond")


class MolfileParserException(Exception):
    pass
