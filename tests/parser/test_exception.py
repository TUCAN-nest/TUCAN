import re
import pytest
from tucan.io import TucanParserException, graph_from_tucan
from tucan.test_utils import TUCAN_VERSION


@pytest.mark.parametrize(
    "tucan, expected_error_msg",
    [
        # lexer error
        (
            f"{TUCAN_VERSION}CXyz/",
            re.escape(
                f"line 1:13 token recognition error at: 'Xy'\n{TUCAN_VERSION}CXyz/\n{' ' * 13}^"
            ),
        ),
        # parser errors
        (
            f"{TUCAN_VERSION}CH4/(1-2)/(1:=14)",
            re.escape(
                f"line 1:25 missing {{'mass', 'rad'}} at '='\n{TUCAN_VERSION}CH4/(1-2)/(1:=14)\n{' ' * 25}^"
            ),
        ),
        (
            f"{TUCAN_VERSION}CH4/(1-2)/(1:massmass=14)",
            re.escape(
                "line 1:29 extraneous input 'mass' expecting"
                f" '='\n{TUCAN_VERSION}CH4/(1-2)/(1:massmass=14)\n{' ' * 29}^^^^"
            ),
        ),
        (
            f"{TUCAN_VERSION}CH4",
            f"line 1:15 missing '/' at '<EOF>'\n{TUCAN_VERSION}CH4",
        ),
        (
            f"{TUCAN_VERSION}CH4/(1-)",
            re.escape(
                "line 1:19 mismatched input ')' expecting {'1', '2', '3', '4', '5', '6',"
                f" '7', '8', '9', GREATER_THAN_NINE}}\n{TUCAN_VERSION}CH4/(1-)\n{' ' * 19}^"
            ),
        ),
    ],
)
def test_graph_from_tucan_error_msg(tucan, expected_error_msg):
    with pytest.raises(TucanParserException, match=expected_error_msg):
        graph_from_tucan(tucan)
