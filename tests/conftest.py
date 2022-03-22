"""Run tests from root of repository with `python -m pytest -v`"""

import pytest
from pathlib import Path
from networkx.algorithms.components import is_connected
from tucan.io import graph_from_file


def pytest_addoption(parser):
    parser.addoption(
        "--testset", action="store", default="mol", help="testset: one of {mol, cfi}"
    )
    parser.addoption(
        "--only-connected-graphs", action="store_true"
    )  # defaults to False


def pytest_configure(config):
    testset_cmd = config.getoption("--testset")
    only_connected_graphs_cmd = config.getoption("--only-connected-graphs")
    testset = None
    if testset_cmd == "mol":
        testset = list(Path("tests/molfiles").glob("*/*.mol"))
    elif testset_cmd == "cfi":
        testset = list(Path("tests/cfi_rigid_benchmark_graphs").glob("*.col"))
    if only_connected_graphs_cmd:
        testset = [m for m in testset if is_connected(graph_from_file(m))]
    pytest.testset = testset

    class Plugin:
        """https://github.com/pytest-dev/pytest/issues/5027"""

        def idtestfile(m):
            """Generate a test ID."""
            return m.stem

        @pytest.fixture(
            params=testset, ids=idtestfile
        )  # automatically runs the test(s) using this fixture on all molecules in `params`
        def m(self, request):
            return graph_from_file(request.param)

    config.pluginmanager.register(Plugin())
