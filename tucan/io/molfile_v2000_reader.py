from tucan.element_properties import ELEMENT_PROPS, detect_hydrogen_isotopes
from tucan.io.exception import MolfileParserException


def graph_props_from_molfile_v2000(lines: list[str]) -> tuple[dict, dict]:
    # counts line: aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
    atom_count = _to_int(lines[3][0:3])  # aaa
    bond_count = _to_int(lines[3][3:6])  # bbb
    atom_lists_count = _to_int(lines[3][6:9])  # lll

    atom_block_offset = 4
    bond_block_offset = atom_block_offset + atom_count
    properties_block_offset = bond_block_offset + atom_lists_count

    atom_props = _parse_atom_block(
        lines[atom_block_offset : atom_block_offset + atom_count]
    )
    bond_props = _parse_bond_block(
        lines[bond_block_offset : bond_block_offset + bond_count], atom_props
    )
    _parse_properties_block(lines[properties_block_offset:], atom_props)

    return atom_props, bond_props


def _parse_atom_block(lines: list[str]) -> dict:
    return {atom_index: _parse_atom_line(line) for atom_index, line in enumerate(lines)}


def _parse_atom_line(line: str) -> dict:
    # atom line: xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
    element_symbol = line[31:34].strip(" ")  # aaa

    element_symbol, isotope_mass = detect_hydrogen_isotopes(element_symbol)

    atom_props = {
        "element_symbol": element_symbol,
        "atomic_number": ELEMENT_PROPS[element_symbol]["atomic_number"],
        "partition": 0,
        "x_coord": _to_float(line[0:10]),  # xxxxx.xxxx
        "y_coord": _to_float(line[10:20]),  # yyyyy.yyyy
        "z_coord": _to_float(line[20:30]),  # zzzzz.zzzz
    }
    atom_props |= _parse_atom_block_charge(line[36:39])  # ccc

    # Field "dd" (mass difference) is ignored. Only consider hydrogen
    # isotopes (D and T) here and "M  ISO" in the properties block (later).
    if isotope_mass:
        atom_props["mass"] = isotope_mass

    return atom_props


def _parse_atom_block_charge(s: str) -> dict:
    match _to_int(s):
        case 1:
            return {"chg": 3}
        case 2:
            return {"chg": 2}
        case 3:
            return {"chg": 1}
        case 4:
            return {"rad": 2}  # doublet radical
        case 5:
            return {"chg": -1}
        case 6:
            return {"chg": -2}
        case 7:
            return {"chg": -3}
        case _:
            return {}  # ignore silently


def _parse_bond_block(lines: list[str], atom_props: dict) -> dict:
    return dict([_parse_bond_line(line, atom_props) for line in lines])


def _parse_bond_line(
    line: str, atom_props: dict
) -> tuple[tuple[int, int], dict[str, int]]:
    # bond line: 111222tttsssxxxrrrccc
    index1 = _to_int(line[0:3]) - 1  # 111
    index2 = _to_int(line[3:6]) - 1  # 222

    _validate_atom_index(index1, atom_props, line)
    _validate_atom_index(index2, atom_props, line)

    bond_props = {"bond_type": _to_int(line[6:9])}  # ttt

    return (index1, index2), bond_props


def _parse_properties_block(lines: list[str], atom_props: dict):
    reset_chg_and_rad = False
    reset_mass = False

    additional_props: dict = {}
    for line in lines:
        if line.startswith("M  CHG"):
            # M  CHGnn8 aaa vvv ...
            _merge_tuples_into_additional_props(
                _parse_atom_value_assignments(line, atom_props), "chg", additional_props
            )
            reset_chg_and_rad = True
        elif line.startswith("M  RAD"):
            # M  RADnn8 aaa vvv ...
            _merge_tuples_into_additional_props(
                _parse_atom_value_assignments(line, atom_props), "rad", additional_props
            )
            reset_chg_and_rad = True
        elif line.startswith("M  ISO"):
            # M  ISOnn8 aaa vvv ...
            _merge_tuples_into_additional_props(
                _parse_atom_value_assignments(line, atom_props),
                "mass",
                additional_props,
            )
            reset_mass = True
        elif line == "M  END":
            break  # else of this for loop is not entered
    else:
        raise MolfileParserException(
            'Could not find end of properties block ("M  END")'
        )

    if reset_chg_and_rad:
        # CHG or RAD lines supersede all charge and radical values from the atom block.
        _clear_atom_prop("chg", atom_props)
        _clear_atom_prop("rad", atom_props)
    if reset_mass:
        # ISO lines supersede all isotope values from the atom block.
        _clear_atom_prop("mass", atom_props)

    _merge_atom_props_and_additional_props(atom_props, additional_props)


def _merge_tuples_into_additional_props(
    atom_index_and_value_tuples: list[tuple[int, int]], key: str, additional_props: dict
):
    for atom_index, value in atom_index_and_value_tuples:
        if atom_index in additional_props:
            additional_props[atom_index][key] = value
        else:
            additional_props[atom_index] = {key: value}


def _parse_atom_value_assignments(line: str, atom_props: dict) -> list[tuple[int, int]]:
    # properties line example: M  CHGnn8 aaa vvv ...
    number_of_entries = _to_int(line[6:9])  # nn8
    tuple_offset = 10
    tuple_length = 8

    assignments = []
    for i in range(number_of_entries):
        tuple_start = tuple_offset + i * tuple_length
        atom_index = _to_int(line[tuple_start : tuple_start + 3]) - 1  # aaa
        value = _to_int(line[tuple_start + 4 : tuple_start + 7])  # vvv

        _validate_atom_index(atom_index, atom_props, line)
        assignments.append((atom_index, value))

    return assignments


def _validate_atom_index(index: int, atom_props: dict, line: str):
    if index not in atom_props:
        raise MolfileParserException(f'Unknown atom index {index + 1} in line "{line}"')


def _clear_atom_prop(key: str, atom_props: dict):
    for atom_prop in atom_props.values():
        # Remove key from dict, but don't raise KeyError if it doesn't exist.
        atom_prop.pop(key, None)


def _merge_atom_props_and_additional_props(atom_props: dict, additional_props: dict):
    for atom_index, props in atom_props.items():
        if atom_index in additional_props:
            props |= additional_props[atom_index]


def _to_int(s: str) -> int:
    if not s.strip(" "):
        return 0
    return int(s)


def _to_float(s: str) -> float:
    if not s.strip(" "):
        return 0
    return float(s)
