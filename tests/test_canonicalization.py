"""Run tests from root of repository with `python -m pytest -v`"""

from canonymous.canonicalization import graph_from_molfile, canonicalize_molecule
from canonymous.utils import permute_molecule, serialize_molecule
from pathlib import Path
import networkx as nx
import random
import pytest


def idtestfile(testfile):
    """Generate a test ID."""
    return testfile.stem

@pytest.fixture(params=list(Path("tests/testfiles").glob("*/*.mol")),
                ids=idtestfile)    # automatically runs the test(s) using this fixture on all values of `params`
def testfiles(request):
    return request.param

def test_permutation(testfiles):
    m = graph_from_molfile(testfiles)
    # Enforce permutation for graphs with at least 2 edges that aren't fully connected (i.e., complete).
    if m.number_of_edges() <= 1 or nx.density(m) == 1:
        return
    permutation_seed = .42
    m_permu = permute_molecule(m, random_seed=permutation_seed)
    assert m.edges != m_permu.edges

def test_invariance(testfiles, n_runs=10, random_seed=random.random()):
    """Eindeutigkeit."""
    m = graph_from_molfile(testfiles)
    root_atom = 0
    random.seed(random_seed)
    for _ in range(n_runs):
        permutation_seed = random.random()
        m_permu = permute_molecule(m, random_seed=permutation_seed)
        m_canon = canonicalize_molecule(m, root_atom)
        m_permu_canon = canonicalize_molecule(m_permu, root_atom)
        m_serialized = serialize_molecule(m_canon)
        m_permu_serialized = serialize_molecule(m_permu_canon)
        assert m_serialized == m_permu_serialized

def test_bijection():
    """Eineindeutigkeit."""
    testfiles = list(Path("tests/testfiles").glob("*/*.mol"))
    serializations = set()
    for f in testfiles:
        m = graph_from_molfile(f)
        m_serialized = serialize_molecule(canonicalize_molecule(m, 0))
        assert m_serialized not in serializations, f"duplicate: {f.stem}"
        serializations.add(m_serialized)

def test_root_atom_independence(testfiles):
    """
    Terephthalic acid reveals why exhaustive permutation of partitions
    (CCAP step in Ivaanciuc, https://doi.org/10.1002/9783527618279.ch7a) is
    necessary.
    However, contrary to Ivanciuc's assertion, not all partitions must be
    permuted. Only those partitions must be permuted that have multiple
    neighbors from the same partition (see graph of
    terephthalic acid in conjunction with output of failed test).
    """
    m = graph_from_molfile(testfiles)
    m_canon = canonicalize_molecule(m, 0)
    for root_atom in range(1, m.number_of_nodes()):
        assert m_canon.edges == canonicalize_molecule(m, root_atom).edges
