import pytest
from antlr4.error.Errors import ParseCancellationException
from tucan.element_properties import element_symbols

from tucan.parser.parser import _prepare_parser


def _parse_sum_formula(s):
    parser = _prepare_parser(s)

    # invoke the parser on rule "sum_formula_start"
    parser.sum_formula_start()


@pytest.mark.parametrize(
    "sum_formula",
    [
        "CHCl3",
        "ClH",
        "H2",
        "C123456789Zr987654321",
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
    ],
)
def test_cannot_parse_sum_formula(sum_formula):
    with pytest.raises(ParseCancellationException):
        _parse_sum_formula(sum_formula)


@pytest.mark.parametrize(
    "symbol", element_symbols
)
def test_grammar_has_all_element_symbols(symbol):
    try:
        _parse_sum_formula(symbol)
    except ParseCancellationException:
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
