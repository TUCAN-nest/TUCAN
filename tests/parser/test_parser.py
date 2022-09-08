import pytest

from tucan.parser.parser import _prepare_parser, _walk_tree


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
