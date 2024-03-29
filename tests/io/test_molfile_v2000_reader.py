import pytest
from pathlib import Path

from tucan.graph_attributes import (
    ATOMIC_NUMBER,
    CHG,
    ELEMENT_SYMBOL,
    MASS,
    PARTITION,
    RAD,
    X_COORD,
    Y_COORD,
    Z_COORD,
)
from tucan.io import graph_from_file
from tucan.io.exception import MolfileParserException
from tucan.io.molfile_v2000_reader import (
    _merge_tuples_into_additional_attributes,
    _parse_atom_line,
    _parse_atom_value_assignments,
    _to_int,
    _to_float,
)


@pytest.mark.parametrize(
    "mol",
    list(map(lambda path: path.stem, Path("tests/molfiles_v2000").glob("*/*.mol"))),
)
def test_graphs_from_v2000_and_v3000_molfiles_match(mol):
    graph_v2000 = graph_from_file(f"tests/molfiles_v2000/{mol}/{mol}.mol")
    graph_v3000 = graph_from_file(f"tests/molfiles/{mol}/{mol}.mol")

    assert graph_v2000.nodes(data=True) == graph_v3000.nodes(data=True)

    # NetworkX's EdgeDataView objects cannot be matched directly
    for e1, e2 in zip(
        graph_v2000.edges(data=True), graph_v3000.edges(data=True), strict=True
    ):
        assert e1 == e2


@pytest.mark.parametrize(
    "line, expected_additional_attr",
    [
        (
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
            {},
        ),
        (
            "    0.0000    0.0000    0.0000 C   0  1  0  0  0  0  0  0  0  0  0  0",
            {CHG: 3},
        ),
        (
            "    0.0000    0.0000    0.0000 C   0  2  0  0  0  0  0  0  0  0  0  0",
            {CHG: 2},
        ),
        (
            "    0.0000    0.0000    0.0000 C   0  3  0  0  0  0  0  0  0  0  0  0",
            {CHG: 1},
        ),
        (
            "    0.0000    0.0000    0.0000 C   0  4  0  0  0  0  0  0  0  0  0  0",
            {RAD: 2},
        ),
        (
            "    0.0000    0.0000    0.0000 C   0  5  0  0  0  0  0  0  0  0  0  0",
            {CHG: -1},
        ),
        (
            "    0.0000    0.0000    0.0000 C   0  6  0  0  0  0  0  0  0  0  0  0",
            {CHG: -2},
        ),
        (
            "    0.0000    0.0000    0.0000 C   0  7  0  0  0  0  0  0  0  0  0  0",
            {CHG: -3},
        ),
        (
            "    0.0000    0.0000    0.0000 C   0  8  0  0  0  0  0  0  0  0  0  0",
            {},  # ignored
        ),
    ],
)
def test_parse_atom_line_charge_field(line, expected_additional_attr):
    expected_attrs = {
        ELEMENT_SYMBOL: "C",
        ATOMIC_NUMBER: 6,
        PARTITION: 0,
        X_COORD: 0,
        Y_COORD: 0,
        Z_COORD: 0,
    }
    expected_attrs.update(expected_additional_attr)

    atom_attrs = _parse_atom_line(line)

    assert atom_attrs == expected_attrs


@pytest.mark.parametrize(
    "tuples, additional_attrs, expected_additional_attrs_after_merge",
    [
        (
            [
                (0, 2),
                (1, 13),
                (2, 3),  # atom index is not in additional_attrs yet
            ],
            {
                0: {CHG: 2},  # will add new key
                1: {MASS: 1},  # will overwrite value
            },
            {
                0: {CHG: 2, MASS: 2},
                1: {MASS: 13},
                2: {MASS: 3},
            },
        ),
    ],
)
def test_merge_tuples_into_additional_attributes(
    tuples, additional_attrs, expected_additional_attrs_after_merge
):
    _merge_tuples_into_additional_attributes(tuples, MASS, additional_attrs)
    assert additional_attrs == expected_additional_attrs_after_merge


@pytest.mark.parametrize(
    "line, atom_attrs, expected_tuples",
    [
        (
            "M  CHG  2  10   2   2 -15",
            {
                1: {},
                9: {},
            },
            [
                (9, 2),
                (1, -15),
            ],
        ),
    ],
)
def test_parse_atom_value_assignments(line, atom_attrs, expected_tuples):
    assert _parse_atom_value_assignments(line, atom_attrs) == expected_tuples


@pytest.mark.parametrize(
    "line, atom_attrs, expected_error_msg",
    [
        (
            "M  CHG  1   5   2",
            {},
            r'Unknown atom index 5 in line "M  CHG  1   5   2"',
        ),
    ],
)
def test_parse_atom_value_assignments_unknown_atom_index_raises_exception(
    line, atom_attrs, expected_error_msg
):
    with pytest.raises(
        MolfileParserException,
        match=expected_error_msg,
    ):
        _parse_atom_value_assignments(line, atom_attrs)


@pytest.mark.parametrize(
    "s, expected_result",
    [
        ("  12 ", 12),
        (" -2    ", -2),
        ("   ", 0),
        ("", 0),
    ],
)
def test_to_int(s, expected_result):
    assert _to_int(s) == expected_result


@pytest.mark.parametrize(
    "s, expected_result",
    [
        ("   3.14 ", 3.14),
        (" -123.456    ", -123.456),
        ("   ", 0),
        ("", 0),
    ],
)
def test_to_float(s, expected_result):
    assert _to_float(s) == expected_result
