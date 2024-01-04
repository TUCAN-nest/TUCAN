import pytest

from tucan.graph_attributes import X_COORD, Y_COORD, Z_COORD
from tucan.io import graph_from_molfile_text, graph_to_molfile
from tucan.io.molfile_writer import _add_header, _add_v30_line


@pytest.mark.parametrize(
    "molfile",
    [
        # molfiles with additional atom attributes
        "tests/molfiles/tnt/tnt.mol",  # CHG
        "tests/molfiles/TEMPO/TEMPO.mol",  # RAD
        "tests/molfiles/water-d1_3/water-d1_3.mol",  # MASS
    ],
)
def test_roundtrip_molfile_graph_molfile(molfile):
    with open(molfile) as file:
        molfile = file.read()
        graph = graph_from_molfile_text(molfile)

        # skip Molfile header for assertion
        assert graph_to_molfile(graph).splitlines()[3:] == molfile.splitlines()[3:]


def test_recalculate_atom_coordinates():
    with open("tests/molfiles/tnt/tnt.mol") as file:
        molfile = file.read()
        graph = graph_from_molfile_text(molfile)

        atoms_with_orig_coords = graph_from_molfile_text(graph_to_molfile(graph)).nodes(
            data=True
        )
        atoms_with_new_coords = graph_from_molfile_text(
            graph_to_molfile(graph, calc_coordinates=True)
        ).nodes(data=True)

        for key, orig_atom_attrs in atoms_with_orig_coords:
            assert orig_atom_attrs[X_COORD] != atoms_with_new_coords[key][X_COORD]
            assert orig_atom_attrs[Y_COORD] != atoms_with_new_coords[key][Y_COORD]
            assert atoms_with_new_coords[key][Z_COORD] == 0.0


def test_add_header():
    lines = []
    _add_header(lines)

    assert len(lines) == 4
    assert lines[0] == ""
    assert lines[1][:7] == "  TUCAN"
    assert lines[1][20:22] == "3D"
    assert lines[2] == ""
    assert lines[3] == "  0  0  0     0  0            999 V3000"


@pytest.mark.parametrize(
    "line, expected_lines",
    [
        # line with 24 chars -> no wrapping
        (
            "1 C -2.58261 2.35748 0 0",
            [
                "M  V30 1 C -2.58261 2.35748 0 0",
            ],
        ),
        # line with 72 chars -> no wrapping
        (
            "1" * 72,
            [
                "M  V30 " + "1" * 72,
            ],
        ),
        # line with 73 chars -> wrap
        (
            "1" * 71 + "2" * 2,
            [
                "M  V30 " + "1" * 71 + "-",
                "M  V30 " + "2" * 2,
            ],
        ),
        # wrap into several lines
        (
            "1" * 71 + "2" * 71 + "3" * 71 + "4" * 20,
            [
                "M  V30 " + "1" * 71 + "-",
                "M  V30 " + "2" * 71 + "-",
                "M  V30 " + "3" * 71 + "-",
                "M  V30 " + "4" * 20,
            ],
        ),
    ],
)
def test_add_v30_line(line, expected_lines):
    lines = []
    _add_v30_line(lines, line)
    assert lines == expected_lines
