import inspect
import pytest
from tucan.parser.parser import TucanParserException, parse_tucan


@pytest.mark.parametrize(
    "tucan, expected_error_msg",
    [
        # lexer error
        (
            "CXyz/",
            """line 1:1 token recognition error at: 'Xy'
            CXyz/
             ^""",
        ),
        # parser errors
        (
            "CH4/(1-2)/(1:=14)",
            """line 1:13 missing {'mass', 'rad'} at '='
            CH4/(1-2)/(1:=14)
                         ^""",
        ),
        (
            "CH4/(1-2)/(1:massmass=14)",
            """line 1:17 extraneous input 'mass' expecting '='
            CH4/(1-2)/(1:massmass=14)
                             ^^^^""",
        ),
        (
            "CH4",
            """line 1:3 missing '/' at '<EOF>'
            CH4
               """,
        ),
        (
            "CH4/(1-)",
            """line 1:7 mismatched input ')' expecting {'1', '2', '3', '4', '5', '6', '7', '8', '9', GREATER_THAN_NINE}
            CH4/(1-)
                   ^""",
        ),
    ],
)
def test_parse_tucan_error_msg(tucan, expected_error_msg):
    with pytest.raises(TucanParserException) as excinfo:
        parse_tucan(tucan)
    assert str(excinfo.value) == inspect.cleandoc(expected_error_msg)
