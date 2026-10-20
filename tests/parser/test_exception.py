import re
import pytest
from tucan.io import TucanParserException, graph_from_tucan
from tucan.test_utils import TUCAN_VERSION


@pytest.mark.parametrize(
    "tucan, expected_error_msg",
    [
        (
            f"{TUCAN_VERSION}CXyz/",
            re.escape(
                f"line 1:13 unknown element symbol at 'X'\n{TUCAN_VERSION}CXyz/\n{' ' * 13}^"
            ),
        ),
        (
            f"{TUCAN_VERSION}CH4/(1-2)/(1:=14)",
            re.escape(
                f"line 1:25 expected 'mass' or 'rad' but found '='\n{TUCAN_VERSION}CH4/(1-2)/(1:=14)\n{' ' * 25}^"
            ),
        ),
        (
            f"{TUCAN_VERSION}CH4/(1-2)/(1:massmass=14)",
            re.escape(
                f"line 1:29 expected '=' but found 'm'\n{TUCAN_VERSION}CH4/(1-2)/(1:massmass=14)\n{' ' * 29}^"
            ),
        ),
        (
            f"{TUCAN_VERSION}CH4",
            re.escape(
                f"line 1:15 expected '/' but found <EOF>\n{TUCAN_VERSION}CH4\n{' ' * 15}^"
            ),
        ),
        (
            f"{TUCAN_VERSION}CH4/(1-)",
            re.escape(
                f"line 1:19 expected a number but found ')'\n{TUCAN_VERSION}CH4/(1-)\n{' ' * 19}^"
            ),
        ),
    ],
)
def test_graph_from_tucan_error_msg(tucan, expected_error_msg):
    with pytest.raises(TucanParserException, match=expected_error_msg):
        graph_from_tucan(tucan)
