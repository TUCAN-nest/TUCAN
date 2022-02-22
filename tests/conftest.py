import pytest
from pathlib import Path


def pytest_addoption(parser):
    parser.addoption(
        "--testset", action="store", default="mol", help="testset: one of {mol, cfi}"
    )


def pytest_configure(config):
    testset_cmd = config.getoption("--testset")
    pytest.testset = None
    if testset_cmd == "mol":
        pytest.testset = list(Path("tests/molfiles").glob("*/*.mol"))
    elif testset_cmd == "cfi":
        pytest.testset = list(Path("tests/cfi_rigid_benchmark_graphs").glob("*.col"))
