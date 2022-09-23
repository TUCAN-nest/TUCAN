"""Run tests from root of repository with `pytest -v`."""

import pytest
import networkx as nx
from pathlib import Path
from networkx.algorithms.components import is_connected
from tucan.io import graph_from_file


def graph_from_dimacs(filepath):
    with open(filepath) as file:
        filecontent = file.read()
    lines = [line.rstrip().split(" ") for line in filecontent.splitlines()]
    lines = [line for line in lines if line[0] in ["p", "e"]]
    atom_count = int(lines[0][2])
    node_labels = range(atom_count)
    bonds = [
        (int(line[1]) - 1, int(line[2]) - 1) for line in lines[1:]
    ]  # make bond-indices zero-based

    graph = nx.Graph()
    graph.add_nodes_from(node_labels)
    nx.set_node_attributes(graph, "C", "element_symbol")
    nx.set_node_attributes(graph, 1, "atomic_number")
    graph.add_edges_from(bonds)

    return graph


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
    filereader = None
    if testset_cmd == "mol":
        testset = list(Path("tests/molfiles").glob("*/*.mol"))
        filereader = graph_from_file
    elif testset_cmd == "cfi":
        testset = list(Path("tests/cfi_rigid_benchmark_graphs").glob("*.col"))
        filereader = graph_from_dimacs
    if only_connected_graphs_cmd:
        testset = [m for m in testset if is_connected(filereader(m))]
    pytest.testset = testset
    pytest.filereader = filereader

    class Plugin:
        """https://github.com/pytest-dev/pytest/issues/5027"""

        def idtestfile(m):
            """Generate a test ID."""
            return m.stem

        @pytest.fixture(
            params=testset, ids=idtestfile
        )  # automatically runs the test(s) using this fixture on all molecules in `params`
        def m(self, request):
            # see https://stackoverflow.com/a/28198398
            if skip_ids_marker := request.node.get_closest_marker("skip_ids"):
                test_id = request.node.callspec.id
                if test_id in skip_ids_marker.args[0]:
                    pytest.skip()

            return filereader(request.param)

    config.pluginmanager.register(Plugin())
    config.addinivalue_line(
        "markers", "skip_ids: ignore tests that match the given ids"
    )
