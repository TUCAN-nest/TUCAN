import re
import pytest
from tucan.io import TucanParserException, graph_from_tucan


@pytest.mark.parametrize(
    "tucan, expected_error_msg",
    [
        # lexer error
        (
            "CXyz/",
            re.escape("line 1:1 token recognition error at: 'Xy'\nCXyz/\n ^"),
        ),
        # parser errors
        (
            "CH4/(1-2)/(1:=14)",
            re.escape(
                "line 1:13 missing {'mass', 'rad'} at '='\nCH4/(1-2)/(1:=14)\n             ^"
            ),
        ),
        (
            "CH4/(1-2)/(1:massmass=14)",
            re.escape(
                "line 1:17 extraneous input 'mass' expecting '='\nCH4/(1-2)/(1:massmass=14)\n                 ^^^^"
            ),
        ),
        (
            "CH4",
            "line 1:3 missing '/' at '<EOF>'\nCH4",
        ),
        (
            "CH4/(1-)",
            re.escape(
                "line 1:7 mismatched input ')' expecting {'1', '2', '3', '4', '5', '6', '7', '8', '9', GREATER_THAN_NINE}\nCH4/(1-)\n       ^"
            ),
        ),
    ],
)
def test_graph_from_tucan_error_msg(tucan, expected_error_msg):
    with pytest.raises(TucanParserException, match=expected_error_msg):
        graph_from_tucan(tucan)
