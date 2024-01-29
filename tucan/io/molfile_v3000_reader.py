import re
from collections import deque
from typing import Any
from tucan.element_attributes import ELEMENT_ATTRS, detect_hydrogen_isotopes
from tucan.graph_attributes import (
    ATOMIC_NUMBER,
    BOND_TYPE,
    CHG,
    ELEMENT_SYMBOL,
    MASS,
    PARTITION,
    RAD,
    X_COORD,
    Y_COORD,
    Z_COORD,
)
from tucan.io.exception import MolfileParserException


def graph_attributes_from_molfile_v3000(
    lines: list[str],
) -> tuple[dict[int, dict[str, Any]], dict[tuple[int, int], dict[str, int]]]:
    tokenized_lines = _tokenize_lines(lines)

    _validate_counts_line(tokenized_lines)
    atom_attrs, star_atoms = _parse_atom_block(tokenized_lines)
    bond_attrs = _parse_bond_block(tokenized_lines, star_atoms)
    _validate_bond_indices(bond_attrs, atom_attrs)

    return atom_attrs, bond_attrs


def _tokenize_lines(lines: list[str]) -> list[list[str]]:
    lines = _concat_lines_with_dash(lines)
    split_lines = [line.rstrip().split(" ") for line in lines]

    return [[value for value in line if value != ""] for line in split_lines]


def _concat_lines_with_dash(line: list[str]) -> list[str]:
    final_lines = []
    line_deque = deque(line)

    while line_deque:
        curr_line = line_deque.popleft()

        if not line_deque:
            final_lines.append(curr_line)
            break

        if not (curr_line.startswith("M  V30 ") and curr_line.endswith("-")):
            final_lines.append(curr_line)
            continue

        next_line = line_deque.popleft()
        if not next_line.startswith("M  V30 "):
            raise MolfileParserException(
                f'Invalid concatenation of lines "{curr_line}" and "{next_line}"'
            )

        line_deque.appendleft(curr_line[0:-1] + next_line[7:])

    return final_lines


def _validate_counts_line(lines: list[list[str]]) -> None:
    if lines[5][2] != "COUNTS" or len(lines[5]) < 5:
        badline = " ".join(lines[5])

        raise MolfileParserException(f'Bad counts line: "{badline}"')


def _parse_atom_block(
    lines: list[list[str]],
) -> tuple[dict[int, dict[str, Any]], list[int]]:
    atom_count = int(lines[5][3])
    atom_block_offset = 7

    if (begin_atom_str := " ".join(lines[atom_block_offset - 1][2:])) != "BEGIN ATOM":
        raise MolfileParserException(
            f'Expected "BEGIN ATOM" on line {atom_block_offset}, found'
            f' "{begin_atom_str}"'
        )
    if (
        end_atom_str := " ".join(lines[atom_block_offset + atom_count][2:])
    ) != "END ATOM":
        raise MolfileParserException(
            f'Expected "END ATOM" on line {atom_block_offset + atom_count + 1}, found'
            f' "{end_atom_str}"'
        )

    atom_attrs = {}
    star_atoms = []
    for line in lines[atom_block_offset : atom_block_offset + atom_count]:
        atom_index = int(line[2]) - 1
        attrs, is_star_atom = _parse_atom_attributes(line)
        if is_star_atom:
            star_atoms.append(atom_index)
            continue
        atom_attrs[atom_index] = attrs

    return atom_attrs, star_atoms


def _parse_atom_attributes(
    line: list[str],
) -> tuple[dict[str, Any], bool]:
    element_symbol = line[3]
    if element_symbol == "*":
        return {}, True

    element_symbol, isotope_mass = detect_hydrogen_isotopes(element_symbol)

    atom_attrs = {
        ELEMENT_SYMBOL: element_symbol,
        ATOMIC_NUMBER: ELEMENT_ATTRS[element_symbol][ATOMIC_NUMBER],
        PARTITION: 0,
        X_COORD: float(line[4]),
        Y_COORD: float(line[5]),
        Z_COORD: float(line[6]),
    }

    optional_attrs = {
        CHG: [int(i.split("=")[1]) for i in line if "CHG" in i],
        MASS: (
            [int(i.split("=")[1]) for i in line if "MASS" in i]
            if not isotope_mass
            else [isotope_mass]
        ),
        RAD: [int(i.split("=")[1]) for i in line if "RAD" in i],
    }
    for key, val in optional_attrs.items():
        if val:
            atom_attrs[key] = val.pop()

    return atom_attrs, False


def _parse_bond_block(
    lines: list[list[str]], star_atoms: list[int]
) -> dict[tuple[int, int], dict[str, int]]:
    atom_count = int(lines[5][3])
    bond_count = int(lines[5][4])

    # bond block is optional
    if bond_count == 0:
        return {}

    atom_block_offset = 7
    bond_block_offset = atom_block_offset + atom_count + 2

    if (begin_bond_str := " ".join(lines[bond_block_offset - 1][2:])) != "BEGIN BOND":
        raise MolfileParserException(
            f'Expected "BEGIN BOND" on line {bond_block_offset}, found'
            f' "{begin_bond_str}"'
        )
    if (
        end_bond_str := " ".join(lines[bond_block_offset + bond_count][2:])
    ) != "END BOND":
        raise MolfileParserException(
            f'Expected "END BOND" on line {bond_block_offset + bond_count + 1}, found'
            f' "{end_bond_str}"'
        )

    bonds = {}
    for line in lines[bond_block_offset : bond_block_offset + bond_count]:
        atom1_index = int(line[4]) - 1
        atom2_index = int(line[5]) - 1
        bond_attrs = _parse_bond_attributes(line)

        atom1_is_star = atom1_index in star_atoms
        atom2_is_star = atom2_index in star_atoms
        if atom1_is_star and atom2_is_star:
            raise MolfileParserException(
                f'Two "star" atoms (index {atom1_index + 1} and {atom2_index + 1}) may'
                " not be connected"
            )
        elif atom1_is_star:
            bond_tuples = _parse_bond_line_with_star_atom(line, atom2_index)
        elif atom2_is_star:
            bond_tuples = _parse_bond_line_with_star_atom(line, atom1_index)
        else:
            bond_tuples = [(atom1_index, atom2_index)]

        for t in bond_tuples:
            bonds[t] = bond_attrs.copy()

    return bonds


def _parse_bond_attributes(line: list[str]) -> dict[str, int]:
    return {BOND_TYPE: int(line[3])}


def _parse_bond_line_with_star_atom(
    line: list[str], start_atom_index: int
) -> list[tuple[int, int]]:
    endpts_pattern = re.compile(r"ENDPTS=\(.+\)")
    endpts_match = endpts_pattern.search(" ".join(line))
    if endpts_match is None:
        # silently ignore everything that has no ENDPTS (e.g. use of star atoms in polymers)
        return []
    endpts_token = endpts_match.group()
    numbers = [int(num) for num in endpts_token[8:-1].split()]
    if (expected_n_endpts := numbers[0]) != (n_endpts := len(numbers) - 1):
        raise MolfileParserException(
            f'Error in "{endpts_token}": Expected {expected_n_endpts} endpoints, found'
            f" {n_endpts}"
        )

    return [(start_atom_index, end_atom_index - 1) for end_atom_index in numbers[1:]]


def _validate_bond_indices(
    bond_attrs: dict[tuple[int, int], dict[str, int]],
    atom_attrs: dict[int, dict[str, Any]],
) -> None:
    for bond in bond_attrs.keys():
        _validate_atom_index(bond[0], atom_attrs)
        _validate_atom_index(bond[1], atom_attrs)


def _validate_atom_index(index: int, atom_attrs: dict[int, dict[str, Any]]) -> None:
    if index not in atom_attrs:
        raise MolfileParserException(f"Unknown atom index {index + 1} in bond")
