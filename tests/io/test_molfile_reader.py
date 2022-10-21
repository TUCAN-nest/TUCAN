import pytest
import re
from tucan.io import graph_from_molfile_text, MolfileParserException
from tucan.io.molfile_reader import (
    _concat_lines_with_dash,
    _parse_atom_block_molfile3000,
    _parse_bond_block_molfile3000,
    _read_file,
    graph_from_file,
)


def test_parsing_atom_block():
    filecontent = _read_file("tests/molfiles/tnt/tnt.mol")
    atom_props, star_atoms = _parse_atom_block_molfile3000(filecontent)
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
    assert len(star_atoms) == 0


def test_graph_from_file_with_multi_attachment():
    graph = graph_from_file(
        "tests/molfiles/chromocene-multi-attachment/chromocene-multi-attachment.mol"
    )
    node_indices, elements = zip(*graph.nodes.data("element_symbol"))

    # "star" atoms do not end up as graph nodes
    assert "".join(elements) == 10 * "C" + "Cr" + 10 * "H"

    # node indices are renamed to achieve sequential order
    assert node_indices == tuple(range(21))

    assert list(graph.edges(data=True)) == [
        (0, 10, {"bond_type": 1}),
        (0, 1, {"bond_type": 4}),
        (0, 4, {"bond_type": 4}),
        (0, 12, {"bond_type": 1}),
        (1, 10, {"bond_type": 1}),
        (1, 2, {"bond_type": 4}),
        (1, 11, {"bond_type": 1}),
        (2, 10, {"bond_type": 1}),
        (2, 3, {"bond_type": 4}),
        (2, 20, {"bond_type": 1}),
        (3, 10, {"bond_type": 1}),
        (3, 4, {"bond_type": 4}),
        (3, 14, {"bond_type": 1}),
        (4, 10, {"bond_type": 1}),
        (4, 13, {"bond_type": 1}),
        (5, 6, {"bond_type": 4}),
        (5, 9, {"bond_type": 4}),
        (5, 10, {"bond_type": 1}),
        (5, 15, {"bond_type": 1}),
        (6, 7, {"bond_type": 4}),
        (6, 10, {"bond_type": 1}),
        (6, 18, {"bond_type": 1}),
        (7, 8, {"bond_type": 4}),
        (7, 10, {"bond_type": 1}),
        (7, 19, {"bond_type": 1}),
        (8, 9, {"bond_type": 4}),
        (8, 10, {"bond_type": 1}),
        (8, 17, {"bond_type": 1}),
        (9, 10, {"bond_type": 1}),
        (9, 16, {"bond_type": 1}),
    ]


@pytest.mark.parametrize(
    "molfile, expected_error_msg",
    [
        # missing "BEGIN ATOM"
        (
            "\n\n\n  0  0  0     0  0            999 V3000\n"
            "M  V30 BEGIN CTAB\n"
            "M  V30 COUNTS 1 1\n"
            "M  V30 1 H 0 0 0 0\n",
            'Expected "BEGIN ATOM" on line 7, found "1 H 0 0 0 0"',
        ),
        # missing "END ATOM"
        (
            "\n\n\n  0  0  0     0  0            999 V3000\n"
            "M  V30 BEGIN CTAB\n"
            "M  V30 COUNTS 2 0\n"
            "M  V30 BEGIN ATOM\n"
            "M  V30 1 H 0 0 0 0\n"
            "M  V30 2 H 0 0 0 0\n"
            "M  V30 3 H 0 0 0 0\n",
            'Expected "END ATOM" on line 10, found "3 H 0 0 0 0"',
        ),
    ],
)
def test_parse_atom_block_molfile3000_raises_exception(molfile, expected_error_msg):
    with pytest.raises(
        MolfileParserException,
        match=expected_error_msg,
    ):
        graph_from_molfile_text(molfile)


def test_parsing_bond_block():
    filecontent = _read_file("tests/molfiles/tnt/tnt.mol")
    bond_props = _parse_bond_block_molfile3000(filecontent, [])
    assert bond_props == {
        (0, 4): {"bond_type": 1},
        (0, 2): {"bond_type": 2},
        (0, 15): {"bond_type": 1},
        (1, 2): {"bond_type": 1},
        (1, 3): {"bond_type": 2},
        (1, 16): {"bond_type": 1},
        (2, 6): {"bond_type": 1},
        (3, 5): {"bond_type": 1},
        (3, 12): {"bond_type": 1},
        (4, 5): {"bond_type": 2},
        (4, 9): {"bond_type": 1},
        (5, 17): {"bond_type": 1},
        (6, 8): {"bond_type": 2},
        (6, 7): {"bond_type": 1},
        (9, 11): {"bond_type": 2},
        (9, 10): {"bond_type": 1},
        (12, 14): {"bond_type": 2},
        (12, 13): {"bond_type": 1},
    }


@pytest.mark.parametrize(
    "molfile, expected_error_msg",
    [
        # missing "BEGIN BOND"
        (
            "\n\n\n  0  0  0     0  0            999 V3000\n"
            "M  V30 BEGIN CTAB\n"
            "M  V30 COUNTS 2 1\n"
            "M  V30 BEGIN ATOM\n"
            "M  V30 1 H 0 0 0 0\n"
            "M  V30 2 H 0 0 0 0\n"
            "M  V30 END ATOM\n"
            "M  V30 1 1 1 2",
            'Expected "BEGIN BOND" on line 11, found "1 1 1 2"',
        ),
        # missing "END BOND"
        (
            "\n\n\n  0  0  0     0  0            999 V3000\n"
            "M  V30 BEGIN CTAB\n"
            "M  V30 COUNTS 2 2\n"
            "M  V30 BEGIN ATOM\n"
            "M  V30 1 H 0 0 0 0\n"
            "M  V30 2 H 0 0 0 0\n"
            "M  V30 END ATOM\n"
            "M  V30 BEGIN BOND\n"
            "M  V30 1 1 1 2\n"
            "M  V30 1 1 2 1\n"
            "M  V30 1 1 1 1",
            'Expected "END BOND" on line 14, found "1 1 1 1"',
        ),
    ],
)
def test_parse_bond_block_molfile3000_raises_exception(molfile, expected_error_msg):
    with pytest.raises(
        MolfileParserException,
        match=expected_error_msg,
    ):
        graph_from_molfile_text(molfile)


@pytest.mark.parametrize(
    "lines, expected_lines",
    [
        # no dash, no concat
        (
            [
                "M  V30 1 H 0 0 0 0",
            ],
            [
                "M  V30 1 H 0 0 0 0",
            ],
        ),
        (
            [
                "M  V30 1 H 0 0 0 0",
                "M  V30 2 O 0 0 0 0",
                "M  V30 3 H 0 0 0 0",
            ],
            [
                "M  V30 1 H 0 0 0 0",
                "M  V30 2 O 0 0 0 0",
                "M  V30 3 H 0 0 0 0",
            ],
        ),
        # non-CTAB lines are not concatenated
        (
            [
                "abc-",
                "M  V30 def",
            ],
            [
                "abc-",
                "M  V30 def",
            ],
        ),
        # example from BIOVIA's CTFile specification
        (
            [
                'M  V30 10 20 30 "abc-',
                'M  V30 def"',
            ],
            [
                'M  V30 10 20 30 "abcdef"',
            ],
        ),
        # dash in last line is silently ignored
        (
            [
                "M  V30 abc-",
            ],
            [
                "M  V30 abc-",
            ],
        ),
        (
            [
                "M  V30 1 H 0 0 0 0 -",
                "M  V30 MASS=2 -",
                "M  V30 RAD=1 -",
            ],
            [
                "M  V30 1 H 0 0 0 0 MASS=2 RAD=1 -",
            ],
        ),
    ],
)
def test_concat_lines_with_dash(lines, expected_lines):
    assert _concat_lines_with_dash(lines) == expected_lines


@pytest.mark.parametrize(
    "lines, expected_error_msg",
    [
        (
            [
                "M  V30 abc-",
                "def",
            ],
            'Invalid concatenation of lines "M  V30 abc-" and "def"',
        ),
    ],
)
def test_concat_lines_with_dash_raises_exception(lines, expected_error_msg):
    with pytest.raises(MolfileParserException, match=expected_error_msg):
        _concat_lines_with_dash(lines)


@pytest.mark.parametrize(
    "molfile",
    [
        "\n\n\n  0  0  0     0  0            999 V2000",
    ],
)
def test_molfile_with_invalid_version_raises_exception(molfile):
    with pytest.raises(
        MolfileParserException,
        match='Invalid Molfile version: Expected "V3000", found "V2000"',
    ):
        graph_from_molfile_text(molfile)


@pytest.mark.parametrize(
    "molfile, expected_error_msg",
    [
        # missing COUNTS line
        (
            "\n\n\n  0  0  0     0  0            999 V3000\n"
            "M  V30 BEGIN CTAB\n"
            "M  V30 BEGIN ATOM",
            'Bad counts line: "M V30 BEGIN ATOM"',
        ),
        # number of bonds missing
        (
            "\n\n\n  0  0  0     0  0            999 V3000\n"
            "M  V30 BEGIN CTAB\n"
            "M  V30 COUNTS 1",
            'Bad counts line: "M V30 COUNTS 1"',
        ),
    ],
)
def test_molfile_with_invalid_counts_line_raises_exception(molfile, expected_error_msg):
    with pytest.raises(
        MolfileParserException,
        match=expected_error_msg,
    ):
        graph_from_molfile_text(molfile)


@pytest.mark.parametrize(
    "molfile, expected_error_msg",
    [
        (
            "\n\n\n  0  0  0     0  0            999 V3000\n"
            "M  V30 BEGIN CTAB\n"
            "M  V30 COUNTS 2 1\n"
            "M  V30 BEGIN ATOM\n"
            "M  V30 1 H 0 0 0 0\n"
            "M  V30 2 H 0 0 0 0\n"
            "M  V30 END ATOM\n"
            "M  V30 BEGIN BOND\n"
            "M  V30 1 1 3 1\n"
            "M  V30 END BOND",
            "Unknown atom index 3 in bond",
        ),
        (
            "\n\n\n  0  0  0     0  0            999 V3000\n"
            "M  V30 BEGIN CTAB\n"
            "M  V30 COUNTS 2 1\n"
            "M  V30 BEGIN ATOM\n"
            "M  V30 1 H 0 0 0 0\n"
            "M  V30 2 H 0 0 0 0\n"
            "M  V30 END ATOM\n"
            "M  V30 BEGIN BOND\n"
            "M  V30 1 1 2 5\n"
            "M  V30 END BOND",
            "Unknown atom index 5 in bond",
        ),
    ],
)
def test_molfile_with_invalid_bond_index_raises_exception(molfile, expected_error_msg):
    with pytest.raises(
        MolfileParserException,
        match=expected_error_msg,
    ):
        graph_from_molfile_text(molfile)


@pytest.mark.parametrize(
    "molfile, expected_error_msg",
    [
        (
            "\n\n\n  0  0  0     0  0            999 V3000\n"
            "M  V30 BEGIN CTAB\n"
            "M  V30 COUNTS 2 1\n"
            "M  V30 BEGIN ATOM\n"
            "M  V30 1 * 0 0 0 0\n"
            "M  V30 2 * 0 0 0 0\n"
            "M  V30 END ATOM\n"
            "M  V30 BEGIN BOND\n"
            "M  V30 1 1 1 2\n"
            "M  V30 END BOND",
            re.escape('Two "star" atoms (index 1 and 2) may not be connected'),
        ),
        (
            "\n\n\n  0  0  0     0  0            999 V3000\n"
            "M  V30 BEGIN CTAB\n"
            "M  V30 COUNTS 4 1\n"
            "M  V30 BEGIN ATOM\n"
            "M  V30 1 * 0 0 0 0\n"
            "M  V30 2 H 0 0 0 0\n"
            "M  V30 3 H 0 0 0 0\n"
            "M  V30 4 H 0 0 0 0\n"
            "M  V30 END ATOM\n"
            "M  V30 BEGIN BOND\n"
            "M  V30 1 1 1 2 ENDPTS=(1 3 4)\n"
            "M  V30 END BOND",
            re.escape('Error in "ENDPTS=(1 3 4)": Expected 1 endpoints, found 2'),
        ),
    ],
)
def test_molfile_with_star_atoms_raises_exception(molfile, expected_error_msg):
    with pytest.raises(
        MolfileParserException,
        match=expected_error_msg,
    ):
        graph_from_molfile_text(molfile)
