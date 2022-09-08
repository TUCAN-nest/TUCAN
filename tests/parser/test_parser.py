import pytest

from tucan.parser.parser import _prepare_parser, _walk_tree, TucanParserException


def _extract_atoms_from_sum_formula(s):
    parser = _prepare_parser(s)
    tree = parser.sum_formula_start()
    listener = _walk_tree(tree)
    return listener._atoms


@pytest.mark.parametrize(
    "sum_formula, expected_atoms",
    [
        ("", []),
        (
            "CHCl3",
            [
                {"element_symbol": "C", "atomic_number": 6, "partition": 0},
                {"element_symbol": "H", "atomic_number": 1, "partition": 0},
                {"element_symbol": "Cl", "atomic_number": 17, "partition": 0},
                {"element_symbol": "Cl", "atomic_number": 17, "partition": 0},
                {"element_symbol": "Cl", "atomic_number": 17, "partition": 0},
            ],
        ),
        (
            "ClH",
            [
                {"element_symbol": "Cl", "atomic_number": 17, "partition": 0},
                {"element_symbol": "H", "atomic_number": 1, "partition": 0},
            ],
        ),
        (
            "Cu",
            [
                {"atomic_number": 29, "element_symbol": "Cu", "partition": 0},
            ],
        ),
        (
            "CU",
            [
                {"atomic_number": 6, "element_symbol": "C", "partition": 0},
                {"atomic_number": 92, "element_symbol": "U", "partition": 0},
            ],
        ),
    ],
)
def test_atoms_from_sum_formula(sum_formula, expected_atoms):
    assert _extract_atoms_from_sum_formula(sum_formula) == expected_atoms


def _extract_bonds_from_tuples(s):
    parser = _prepare_parser(s)
    tree = parser.tuples_start()
    listener = _walk_tree(tree)
    return listener._bonds


@pytest.mark.parametrize(
    "tuples, expected_bonds",
    [
        ("", []),
        (
            "(1-2)(3-4)(5-6)(7-8)(9-10)(11-12)",
            [(0, 1), (2, 3), (4, 5), (6, 7), (8, 9), (10, 11)],
        ),
        ("(2-1)(1-2)", [(1, 0), (0, 1)]),
        ("(1-2)(1-2)", [(0, 1), (0, 1)]),
        ("(987654321-123456789)", [(987654320, 123456788)]),
    ],
)
def test_bonds_from_tuples(tuples, expected_bonds):
    assert _extract_bonds_from_tuples(tuples) == expected_bonds


def test_self_loop_in_tuples_raises_exception():
    with pytest.raises(TucanParserException) as excinfo:
        _extract_bonds_from_tuples("(1-2)(1-1)(2-2)")
    assert str(excinfo.value) == 'Error in tuple "(1-1)": Self-loops are not allowed.'
