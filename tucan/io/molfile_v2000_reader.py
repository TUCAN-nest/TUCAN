from typing import Any
from tucan.element_attributes import (
    ELEMENT_ATTRS,
    MOLFILE_V2000_CHARGES,
    detect_hydrogen_isotopes,
)
from tucan.graph_attributes import ATOMIC_NUMBER
from tucan.io.exception import MolfileParserException


def graph_attributes_from_molfile_v2000(
    lines: list[str],
) -> tuple[dict[int, dict[str, Any]], dict[tuple[int, int], dict[str, int]]]:
    # counts line: aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
    atom_count = _to_int(lines[3][0:3])  # aaa
    bond_count = _to_int(lines[3][3:6])  # bbb
    atom_lists_count = _to_int(lines[3][6:9])  # lll

    atom_block_offset = 4
    bond_block_offset = atom_block_offset + atom_count
    attribute_block_offset = bond_block_offset + atom_lists_count

    atom_attrs = _parse_atom_block(
        lines[atom_block_offset : atom_block_offset + atom_count]
    )
    bond_attrs = _parse_bond_block(
        lines[bond_block_offset : bond_block_offset + bond_count], atom_attrs
    )
    _parse_attribute_block(lines[attribute_block_offset:], atom_attrs)

    return atom_attrs, bond_attrs


def _parse_atom_block(lines: list[str]) -> dict[int, dict[str, Any]]:
    return {atom_index: _parse_atom_line(line) for atom_index, line in enumerate(lines)}


def _parse_atom_line(line: str) -> dict[str, Any]:
    # atom line: xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
    element_symbol = line[31:34].strip(" ")  # aaa

    element_symbol, isotope_mass = detect_hydrogen_isotopes(element_symbol)

    atom_attrs = {
        "element_symbol": element_symbol,
        ATOMIC_NUMBER: ELEMENT_ATTRS[element_symbol][ATOMIC_NUMBER],
        "partition": 0,
        "x_coord": _to_float(line[0:10]),  # xxxxx.xxxx
        "y_coord": _to_float(line[10:20]),  # yyyyy.yyyy
        "z_coord": _to_float(line[20:30]),  # zzzzz.zzzz
    }
    atom_attrs |= MOLFILE_V2000_CHARGES.get(_to_int(line[36:39]), {})  # ccc

    # Field "dd" (mass difference) is ignored. Only consider hydrogen
    # isotopes (D and T) here and "M  ISO" in the attribute block (later).
    if isotope_mass:
        atom_attrs["mass"] = isotope_mass

    return atom_attrs


def _parse_bond_block(
    lines: list[str], atom_attrs: dict[int, dict[str, Any]]
) -> dict[tuple[int, int], dict[str, int]]:
    return dict([_parse_bond_line(line, atom_attrs) for line in lines])


def _parse_bond_line(
    line: str, atom_attrs: dict[int, dict[str, Any]]
) -> tuple[tuple[int, int], dict[str, int]]:
    # bond line: 111222tttsssxxxrrrccc
    index1 = _to_int(line[0:3]) - 1  # 111
    index2 = _to_int(line[3:6]) - 1  # 222

    _validate_atom_index(index1, atom_attrs, line)
    _validate_atom_index(index2, atom_attrs, line)

    bond_attrs = {"bond_type": _to_int(line[6:9])}  # ttt

    return (index1, index2), bond_attrs


def _parse_attribute_block(
    lines: list[str], atom_attrs: dict[int, dict[str, Any]]
) -> None:
    reset_chg_and_rad = False
    reset_mass = False

    additional_attrs: dict = {}
    for line in lines:
        if line.startswith("M  CHG"):
            # M  CHGnn8 aaa vvv ...
            _merge_tuples_into_additional_attributes(
                _parse_atom_value_assignments(line, atom_attrs), "chg", additional_attrs
            )
            reset_chg_and_rad = True
        elif line.startswith("M  RAD"):
            # M  RADnn8 aaa vvv ...
            _merge_tuples_into_additional_attributes(
                _parse_atom_value_assignments(line, atom_attrs), "rad", additional_attrs
            )
            reset_chg_and_rad = True
        elif line.startswith("M  ISO"):
            # M  ISOnn8 aaa vvv ...
            _merge_tuples_into_additional_attributes(
                _parse_atom_value_assignments(line, atom_attrs),
                "mass",
                additional_attrs,
            )
            reset_mass = True
        elif line == "M  END":
            break  # else of this for loop is not entered
    else:
        raise MolfileParserException('Could not find end of attribute block ("M  END")')

    if reset_chg_and_rad:
        # CHG or RAD lines supersede all charge and radical values from the atom block.
        _clear_atom_attribute("chg", atom_attrs)
        _clear_atom_attribute("rad", atom_attrs)
    if reset_mass:
        # ISO lines supersede all isotope values from the atom block.
        _clear_atom_attribute("mass", atom_attrs)

    _merge_atom_attributes_and_additional_attributes(atom_attrs, additional_attrs)


def _merge_tuples_into_additional_attributes(
    atom_index_and_value_tuples: list[tuple[int, int]], key: str, additional_attrs: dict
) -> None:
    for atom_index, value in atom_index_and_value_tuples:
        if atom_index in additional_attrs:
            additional_attrs[atom_index][key] = value
        else:
            additional_attrs[atom_index] = {key: value}


def _parse_atom_value_assignments(
    line: str, atom_attrs: dict[int, dict[str, Any]]
) -> list[tuple[int, int]]:
    # attribute line example: M  CHGnn8 aaa vvv ...
    number_of_entries = _to_int(line[6:9])  # nn8
    tuple_offset = 10
    tuple_length = 8

    assignments = []
    for i in range(number_of_entries):
        tuple_start = tuple_offset + i * tuple_length
        atom_index = _to_int(line[tuple_start : tuple_start + 3]) - 1  # aaa
        value = _to_int(line[tuple_start + 4 : tuple_start + 7])  # vvv

        _validate_atom_index(atom_index, atom_attrs, line)
        assignments.append((atom_index, value))

    return assignments


def _validate_atom_index(
    index: int, atom_attrs: dict[int, dict[str, Any]], line: str
) -> None:
    if index not in atom_attrs:
        raise MolfileParserException(f'Unknown atom index {index + 1} in line "{line}"')


def _clear_atom_attribute(key: str, atom_attrs: dict[int, dict[str, Any]]) -> None:
    for atom_attr in atom_attrs.values():
        # Remove key from dict, but don't raise KeyError if it doesn't exist.
        atom_attr.pop(key, None)


def _merge_atom_attributes_and_additional_attributes(
    atom_attrs: dict[int, dict[str, Any]], additional_attrs: dict
) -> None:
    for atom_index, attrs in atom_attrs.items():
        if atom_index in additional_attrs:
            attrs |= additional_attrs[atom_index]


def _to_int(s: str) -> int:
    if not s.strip(" "):
        return 0

    return int(s)


def _to_float(s: str) -> float:
    if not s.strip(" "):
        return 0

    return float(s)
