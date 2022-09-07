import pytest
from tucan.canonicalization import canonicalize_molecule
from tucan.element_properties import element_symbols
from tucan.parser.parser import _prepare_parser, TucanParserException, parse_tucan
from tucan.serialization import serialize_molecule


def _parse_sum_formula(s):
    parser = _prepare_parser(s)

    # invoke the parser on rule "sum_formula_start"
    parser.sum_formula_start()


@pytest.mark.parametrize(
    "sum_formula",
    [
        "",
        "CHCl3",
        "ClH",
        "H2",
        "C123456789Zr987654321",
        "Cu",
        "CU",
        "Os",
        "OS",
        "Bi",
        "BI",
        "IIn",
    ],
)
def test_can_parse_sum_formula(sum_formula):
    _parse_sum_formula(sum_formula)


@pytest.mark.parametrize(
    "sum_formula",
    [
        "HCl",
        "H0",
        "H1",
        "H0123",
        "HH",
        "Ar3Ar10",
        "D2O",
        "OT2",
        "h",
        "hCl",
        "Hga",
    ],
)
def test_cannot_parse_sum_formula(sum_formula):
    with pytest.raises(TucanParserException):
        _parse_sum_formula(sum_formula)


@pytest.mark.parametrize("symbol", element_symbols)
def test_grammar_has_all_element_symbols(symbol):
    try:
        _parse_sum_formula(symbol)
    except TucanParserException:
        pytest.fail('Unknown element symbol "' + symbol + '"')


def test_element_order_is_correct():
    elements_without_carbon = element_symbols.copy()
    elements_without_carbon.remove("C")
    elements_without_carbon.sort()
    _parse_sum_formula("".join(elements_without_carbon))

    elements_with_carbon = elements_without_carbon.copy()
    elements_with_carbon.remove("H")
    elements_with_carbon.insert(0, "H")
    elements_with_carbon.insert(0, "C")
    _parse_sum_formula("".join(elements_with_carbon))


def _parse_tuples(s):
    parser = _prepare_parser(s)

    # invoke the parser on rule "tuples_start"
    parser.tuples_start()


@pytest.mark.parametrize(
    "tuples",
    [
        "",
        "(1-2)(3-4)(5-6)(7-8)(9-10)(11-12)",
        "(2-1)(1-2)",
        "(1-2)(1-2)",
        "(1-1)",
        "(987654321-123456789)",
    ],
)
def test_can_parse_tuples(tuples):
    _parse_tuples(tuples)


@pytest.mark.parametrize(
    "tuples",
    [
        "(0-1)",
        "(1-0)",
        "(01-2)",
        "(1-02)",
        "1-2",
        "(1-2",
        "1-2)",
        "(1 2)",
        "(1--2)",
        "(12)",
        "()",
        "(",
        ")",
        "()(1-2)",
        # different dash symbols ;)
        "(1‐2)",
        "(1–2)",
        "(1—2)",
        "(1―2)",
        "(1−2)",
    ],
)
def test_cannot_parse_tuples(tuples):
    with pytest.raises(TucanParserException):
        _parse_tuples(tuples)


def _parse_node_attributes(s):
    parser = _prepare_parser(s)

    # invoke the parser on rule "node_attributes_start"
    parser.node_attributes_start()


@pytest.mark.parametrize(
    "node_attributes",
    [
        "",
        "(1:mass=2)",
        "(2:rad=5)",
        "(3210:mass=10,rad=5)",
        "(1234:rad=5,mass=10)",
        "(1:mass=10)(2:rad=1)",
        "(2:rad=1)(1:mass=10)",
        "(1:mass=123456789)(1:rad=987654321)",
        "(1:mass=10,rad=5)(2:rad=1)(1:rad=3,mass=12)",
        "(2:mass=5,mass=7)",
    ],
)
def test_can_parse_node_attributes(node_attributes):
    _parse_node_attributes(node_attributes)


@pytest.mark.parametrize(
    "node_attributes",
    [
        "()",
        "(1:)",
        "(1:mass=2,)",
        "(2:mass=2,rad=4,)",
        "(01:mass=2)",
        "(:mass=2)",
        "(4321:MAS=2)",
        "(4321:mass2)",
        "(4321:mass 2)",
        "(1:mass=0)",
        "(1:rad=0123)",
        "(1:mass=)",
        "(1:rad)",
        "(2:mass=10",
        "(3:mass=10)(",
        "(",
        ")",
    ],
)
def test_cannot_parse_node_attributes(node_attributes):
    with pytest.raises(TucanParserException):
        _parse_node_attributes(node_attributes)


@pytest.mark.parametrize(
    "tucan",
    [
        "/",
        "//",
        "C2H6O/(1-7)(2-7)(3-7)(4-8)(5-8)(6-9)(7-8)(8-9)",
        "Xe/",
        "C2H4O/(1-5)(2-5)(3-5)(4-7)(5-6)(6-7)/(4:mass=2)(5:mass=14)(6:rad=3)(7:mass=17)",
    ],
)
def test_can_parse_tucan(tucan):
    parse_tucan(tucan)


@pytest.mark.parametrize(
    "tucan",
    [
        "",
        "C2H6O(1-7)(2-7)(3-7)(4-8)(5-8)(6-9)(7-8)(8-9)",
        "C2H6O(1-7)(2-7)(3-7)(4-8)(5-8)(6-9)(7-8)(8-9)/",
    ],
)
def test_cannot_parse_tucan(tucan):
    with pytest.raises(TucanParserException):
        parse_tucan(tucan)


def test_roundtrip_graph_tucan_graph(m):
    m_serialized = serialize_molecule(canonicalize_molecule(m))
    parse_tucan(m_serialized)
