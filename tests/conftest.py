import pytest
from pathlib import Path
from networkx.algorithms.components import is_connected
from tucan.canonicalization import graph_from_file


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
