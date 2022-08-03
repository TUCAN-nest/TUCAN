from typing import List
from tucan.io import (
    _parse_atom_block_molfile3000,
    _parse_bond_block_molfile3000,
    _split_into_tokenized_lines,
)


def test_parsing_atom_block():
    filecontent = _read_file("tests/molfiles/tnt/tnt.mol")
    atom_props = _parse_atom_block_molfile3000(filecontent)
    assert atom_props == {
        0: {
            "element_symbol": "C",
            "atomic_number": 6,
            "partition": 0,
            "x_coord": 11.137,
            "y_coord": -9.481,
            "z_coord": 0.0,
        },
        1: {
            "element_symbol": "C",
            "atomic_number": 6,
            "partition": 0,
            "x_coord": 12.1745,
            "y_coord": -9.4807,
            "z_coord": 0.0,
        },
        2: {
            "element_symbol": "C",
            "atomic_number": 6,
            "partition": 0,
            "x_coord": 11.6567,
            "y_coord": -9.1811,
            "z_coord": 0.0,
        },
        3: {
            "element_symbol": "C",
            "atomic_number": 6,
            "partition": 0,
            "x_coord": 12.1745,
            "y_coord": -10.0809,
            "z_coord": 0.0,
        },
        4: {
            "element_symbol": "C",
            "atomic_number": 6,
            "partition": 0,
            "x_coord": 11.137,
            "y_coord": -10.0835,
            "z_coord": 0.0,
        },
        5: {
            "element_symbol": "C",
            "atomic_number": 6,
            "partition": 0,
            "x_coord": 11.658,
            "y_coord": -10.3804,
            "z_coord": 0.0,
        },
        6: {
            "element_symbol": "N",
            "atomic_number": 7,
            "partition": 0,
            "x_coord": 11.6691,
            "y_coord": -7.3712,
            "z_coord": 0.0,
            "chg": 1,
        },
        7: {
            "element_symbol": "O",
            "atomic_number": 8,
            "partition": 0,
            "x_coord": 12.1887,
            "y_coord": -7.0712,
            "z_coord": 0.0,
            "chg": -1,
        },
        8: {
            "element_symbol": "O",
            "atomic_number": 8,
            "partition": 0,
            "x_coord": 11.1495,
            "y_coord": -7.0712,
            "z_coord": 0.0,
        },
        9: {
            "element_symbol": "N",
            "atomic_number": 7,
            "partition": 0,
            "x_coord": 8.8633,
            "y_coord": -11.1246,
            "z_coord": 0.0,
            "chg": 1,
        },
        10: {
            "element_symbol": "O",
            "atomic_number": 8,
            "partition": 0,
            "x_coord": 9.0299,
            "y_coord": -12.4412,
            "z_coord": 0.0,
            "chg": -1,
        },
        11: {
            "element_symbol": "O",
            "atomic_number": 8,
            "partition": 0,
            "x_coord": 8.3437,
            "y_coord": -10.8246,
            "z_coord": 0.0,
        },
        12: {
            "element_symbol": "N",
            "atomic_number": 7,
            "partition": 0,
            "x_coord": 13.8431,
            "y_coord": -11.1804,
            "z_coord": 0.0,
            "chg": 1,
        },
        13: {
            "element_symbol": "O",
            "atomic_number": 8,
            "partition": 0,
            "x_coord": 14.3627,
            "y_coord": -10.8804,
            "z_coord": 0.0,
            "chg": -1,
        },
        14: {
            "element_symbol": "O",
            "atomic_number": 8,
            "partition": 0,
            "x_coord": 13.3607,
            "y_coord": -12.0324,
            "z_coord": 0.0,
        },
        15: {
            "element_symbol": "H",
            "atomic_number": 1,
            "partition": 0,
            "x_coord": 9.4208,
            "y_coord": -8.4533,
            "z_coord": 0.0,
        },
        16: {
            "element_symbol": "H",
            "atomic_number": 1,
            "partition": 0,
            "x_coord": 14.0661,
            "y_coord": -8.4162,
            "z_coord": 0.0,
        },
        17: {
            "element_symbol": "H",
            "atomic_number": 1,
            "partition": 0,
            "x_coord": 11.2046,
            "y_coord": -12.0581,
            "z_coord": 0.0,
        },
    }


def test_parsing_bond_block():
    filecontent = _read_file("tests/molfiles/tnt/tnt.mol")
    bonds = _parse_bond_block_molfile3000(filecontent)
    assert bonds == [
        (0, 4),
        (2, 0),
        (0, 15),
        (1, 2),
        (3, 1),
        (1, 16),
        (2, 6),
        (5, 3),
        (3, 12),
        (4, 5),
        (4, 9),
        (5, 17),
        (6, 8),
        (6, 7),
        (9, 11),
        (9, 10),
        (12, 14),
        (12, 13),
    ]


def _read_file(filepath: str) -> List[List[str]]:
    with open(filepath) as file:
        filecontent = file.read()
    return _split_into_tokenized_lines(filecontent)