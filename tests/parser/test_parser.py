import pytest
import re

from tucan.graph_attributes import (
    ATOMIC_NUMBER,
    ELEMENT_SYMBOL,
    INVARIANT_CODE,
    MASS,
    PARTITION,
    RAD,
)
from tucan.io import graph_from_tucan, TucanParserException
from tucan.parser.parser import _prepare_parser, _walk_tree
from tucan.test_utils import roundtrip_graph_tucan_graph_tucan_graph


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
                {ELEMENT_SYMBOL: "C", ATOMIC_NUMBER: 6, PARTITION: 0},
                {ELEMENT_SYMBOL: "H", ATOMIC_NUMBER: 1, PARTITION: 0},
                {ELEMENT_SYMBOL: "Cl", ATOMIC_NUMBER: 17, PARTITION: 0},
                {ELEMENT_SYMBOL: "Cl", ATOMIC_NUMBER: 17, PARTITION: 0},
                {ELEMENT_SYMBOL: "Cl", ATOMIC_NUMBER: 17, PARTITION: 0},
            ],
        ),
        (
            "ClH",
            [
                {ELEMENT_SYMBOL: "Cl", ATOMIC_NUMBER: 17, PARTITION: 0},
                {ELEMENT_SYMBOL: "H", ATOMIC_NUMBER: 1, PARTITION: 0},
            ],
        ),
        (
            "Cu",
            [
                {ELEMENT_SYMBOL: "Cu", ATOMIC_NUMBER: 29, PARTITION: 0},
            ],
        ),
        (
            "CU",
            [
                {ELEMENT_SYMBOL: "C", ATOMIC_NUMBER: 6, PARTITION: 0},
                {ELEMENT_SYMBOL: "U", ATOMIC_NUMBER: 92, PARTITION: 0},
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
    expected_error_msg = re.escape(
        'Error in tuple "(1-1)": Self-loops are not allowed.'
    )
    with pytest.raises(TucanParserException, match=expected_error_msg):
        _extract_bonds_from_tuples("(1-2)(1-1)(2-2)")


def _extract_node_attributes(s):
    parser = _prepare_parser(s)
    tree = parser.node_attributes_start()
    listener = _walk_tree(tree)
    return listener._node_attributes


@pytest.mark.parametrize(
    "node_attributes, expected_node_attributes",
    [
        ("", {}),
        ("(1:mass=2)", {0: {MASS: 2}}),
        ("(2:rad=5)", {1: {RAD: 5}}),
        ("(1234:rad=5,mass=10)", {1233: {MASS: 10, RAD: 5}}),
        ("(1:mass=10)(2:rad=1)", {0: {MASS: 10}, 1: {RAD: 1}}),
        (
            "(1:mass=123456789)(1:rad=987654321)",
            {0: {MASS: 123456789, RAD: 987654321}},
        ),
    ],
)
def test_node_attributes(node_attributes, expected_node_attributes):
    assert _extract_node_attributes(node_attributes) == expected_node_attributes


@pytest.mark.parametrize(
    "node_attributes, offending_node_index, offending_key",
    [
        ("(1:mass=10,rad=5)(2:rad=1)(1:rad=3,mass=12)", 1, "rad"),
        ("(2:mass=5,mass=7)", 2, "mass"),
    ],
)
def test_overriding_node_attribute_raises_exception(
    node_attributes, offending_node_index, offending_key
):
    expected_error_msg = (
        f'^Atom {offending_node_index}: Attribute "{offending_key}" was already'
        " defined.$"
    )
    with pytest.raises(TucanParserException, match=expected_error_msg):
        _extract_node_attributes(node_attributes)


@pytest.mark.parametrize(
    "tucan, expected_atoms, expected_bonds",
    [
        ("/", {}, []),
        ("//", {}, []),
        (
            "C2H6O/(1-7)(2-7)(3-7)(4-8)(5-8)(6-9)(7-8)(8-9)",
            {
                0: {
                    ELEMENT_SYMBOL: "H",
                    ATOMIC_NUMBER: 1,
                    PARTITION: 0,
                    INVARIANT_CODE: (1, 0, 0),
                },
                1: {
                    ELEMENT_SYMBOL: "H",
                    ATOMIC_NUMBER: 1,
                    PARTITION: 0,
                    INVARIANT_CODE: (1, 0, 0),
                },
                2: {
                    ELEMENT_SYMBOL: "H",
                    ATOMIC_NUMBER: 1,
                    PARTITION: 0,
                    INVARIANT_CODE: (1, 0, 0),
                },
                3: {
                    ELEMENT_SYMBOL: "H",
                    ATOMIC_NUMBER: 1,
                    PARTITION: 0,
                    INVARIANT_CODE: (1, 0, 0),
                },
                4: {
                    ELEMENT_SYMBOL: "H",
                    ATOMIC_NUMBER: 1,
                    PARTITION: 0,
                    INVARIANT_CODE: (1, 0, 0),
                },
                5: {
                    ELEMENT_SYMBOL: "H",
                    ATOMIC_NUMBER: 1,
                    PARTITION: 0,
                    INVARIANT_CODE: (1, 0, 0),
                },
                6: {
                    ELEMENT_SYMBOL: "C",
                    ATOMIC_NUMBER: 6,
                    PARTITION: 0,
                    INVARIANT_CODE: (6, 0, 0),
                },
                7: {
                    ELEMENT_SYMBOL: "C",
                    ATOMIC_NUMBER: 6,
                    PARTITION: 0,
                    INVARIANT_CODE: (6, 0, 0),
                },
                8: {
                    ELEMENT_SYMBOL: "O",
                    ATOMIC_NUMBER: 8,
                    PARTITION: 0,
                    INVARIANT_CODE: (8, 0, 0),
                },
            },
            [(0, 6), (1, 6), (2, 6), (3, 7), (4, 7), (5, 8), (6, 7), (7, 8)],
        ),
        (
            "Xe/",
            {
                0: {
                    ELEMENT_SYMBOL: "Xe",
                    ATOMIC_NUMBER: 54,
                    PARTITION: 0,
                    INVARIANT_CODE: (54, 0, 0),
                }
            },
            [],
        ),
        (
            "C2H4O/(1-5)(2-5)(3-5)(4-7)(5-6)(6-7)/(4:mass=2)(5:mass=14)(6:rad=3)(7:mass=17)",
            {
                0: {
                    ELEMENT_SYMBOL: "H",
                    ATOMIC_NUMBER: 1,
                    PARTITION: 0,
                    INVARIANT_CODE: (1, 0, 0),
                },
                1: {
                    ELEMENT_SYMBOL: "H",
                    ATOMIC_NUMBER: 1,
                    PARTITION: 0,
                    INVARIANT_CODE: (1, 0, 0),
                },
                2: {
                    ELEMENT_SYMBOL: "H",
                    ATOMIC_NUMBER: 1,
                    PARTITION: 0,
                    INVARIANT_CODE: (1, 0, 0),
                },
                3: {
                    ELEMENT_SYMBOL: "H",
                    ATOMIC_NUMBER: 1,
                    PARTITION: 0,
                    MASS: 2,
                    INVARIANT_CODE: (1, 2, 0),
                },
                4: {
                    ELEMENT_SYMBOL: "C",
                    ATOMIC_NUMBER: 6,
                    PARTITION: 0,
                    MASS: 14,
                    INVARIANT_CODE: (6, 14, 0),
                },
                5: {
                    ELEMENT_SYMBOL: "C",
                    ATOMIC_NUMBER: 6,
                    PARTITION: 0,
                    RAD: 3,
                    INVARIANT_CODE: (6, 0, 3),
                },
                6: {
                    ELEMENT_SYMBOL: "O",
                    ATOMIC_NUMBER: 8,
                    PARTITION: 0,
                    MASS: 17,
                    INVARIANT_CODE: (8, 17, 0),
                },
            },
            [(0, 4), (1, 4), (2, 4), (3, 6), (4, 5), (5, 6)],
        ),
        (
            "CH/(2-1)(1-2)/(2:rad=3,mass=13)",
            {
                0: {
                    ELEMENT_SYMBOL: "H",
                    ATOMIC_NUMBER: 1,
                    PARTITION: 0,
                    INVARIANT_CODE: (1, 0, 0),
                },
                1: {
                    ELEMENT_SYMBOL: "C",
                    ATOMIC_NUMBER: 6,
                    PARTITION: 0,
                    RAD: 3,
                    MASS: 13,
                    INVARIANT_CODE: (6, 13, 3),
                },
            },
            [(0, 1)],
        ),
    ],
)
def test_graph_from_tucan(tucan, expected_atoms, expected_bonds):
    graph = graph_from_tucan(tucan)
    assert dict(graph.nodes(data=True)) == expected_atoms
    assert list(graph.edges) == expected_bonds


def test_roundtrip_molfile_graph_tucan_graph_tucan_graph(m):
    roundtrip_graph_tucan_graph_tucan_graph(m)


@pytest.mark.parametrize(
    "tucan, offending_node_index",
    [
        ("CH/(1-3)", 3),
        ("CH/(1-2)(5-1)", 5),
        ("CH//(3:mass=1)", 3),
        ("CH3//(1:mass=13)(5:rad=3)", 5),
    ],
)
def test_graph_from_tucan_invalid_node_index_raises_exception(
    tucan, offending_node_index
):
    expected_error_msg = f"^Atom with index {offending_node_index} does not exist.$"
    with pytest.raises(TucanParserException, match=expected_error_msg):
        graph = graph_from_tucan(tucan)
