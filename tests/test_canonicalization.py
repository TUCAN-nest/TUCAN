"""Run tests from root of repository with `python -m pytest -v`"""

from canonymous.canonicalization import graph_from_molfile, canonicalize_molecule
from canonymous.utils import permute_molecule
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

def test_invariance(testfiles, n_runs=10, root_atom=None, visualize=False):
    """
    TODO:
    * write separate test for different root atoms (having it in this test clutters test)
    * refactor `permute_molecule`: permutation should be enforced by default (shouldn't be handled in this test)
    * have test emit permutation seed on failure for reproducibility
    """
    m = graph_from_molfile(testfiles)
    # Enforce permutation for graphs with at least 2 edges that aren't fully connected (i.e., complete).
    enforce_permutation = m.number_of_edges() > 1 and nx.density(m) != 1
    root_atoms = [root_atom] * n_runs
    n_nodes = m.number_of_nodes()
    if not root_atom:
        if n_runs <= n_nodes:
            root_atoms = random.sample(range(n_nodes), k=n_runs)    # draw without replacement
        else:
            root_atoms = random.choices(range(n_nodes), k=n_runs)    # draw with replacement

    for i in range(n_runs):
        root_atom = root_atoms[i]
        permutation_seed = random.random()

        m_permu = permute_molecule(m, random_seed=permutation_seed)
        if enforce_permutation:
            while m.edges == m_permu.edges:
                permutation_seed = random.random()
                m_permu = permute_molecule(m, random_seed=permutation_seed)
            assert m.edges != m_permu.edges

        m_canon = canonicalize_molecule(m, root_atom)
        m_permu_canon = canonicalize_molecule(m_permu, root_atom)
        assert m_canon.edges == m_permu_canon.edges
