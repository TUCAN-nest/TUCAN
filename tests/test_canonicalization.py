from tucan.canonicalization import (
    canonicalize_molecule,
    permute_molecule,
    serialize_molecule,
)
from tucan.io import graph_from_file
import networkx as nx
import random
import pytest
from pathlib import Path


def test_permutation(m):
    # Enforce permutation for graphs with at least 2 edges that aren't fully connected (i.e., complete).
    if m.number_of_edges() <= 1:
        pytest.skip("Skipping graph with less than two edges.")
    if nx.density(m) == 1:
        pytest.skip("Skipping fully connected graph.")
    permutation_seed = 0.42
    m_permu = permute_molecule(m, random_seed=permutation_seed)
    assert m.edges != m_permu.edges


def test_invariance(m, n_runs=10, random_seed=random.random(), root_atom=0):
    """Eindeutigkeit."""
    m_canon = canonicalize_molecule(m, root_atom)
    m_serialized = serialize_molecule(m_canon)
    random.seed(random_seed)
    for _ in range(n_runs):
        permutation_seed = random.random()
        m_permu = permute_molecule(m, random_seed=permutation_seed)
        m_permu_canon = canonicalize_molecule(m_permu, root_atom)
        m_permu_serialized = serialize_molecule(m_permu_canon)
        assert m_serialized == m_permu_serialized


def test_bijection():
    """Eineindeutigkeit."""
    serializations = set()
    for f in pytest.testset:
        m = graph_from_file(f)
        m_serialized = serialize_molecule(canonicalize_molecule(m, 0))
        assert m_serialized not in serializations, f"duplicate: {f.stem}"
        serializations.add(m_serialized)


def test_root_atom_independence(m):
    m_canon = canonicalize_molecule(m, 0)
    for root_atom in range(1, m.number_of_nodes()):
        assert m_canon.edges == canonicalize_molecule(m, root_atom).edges


@pytest.mark.parametrize(
    "m, expected_serialization",
    [
        (
            "ferrocene",
            "C10H10Fe/1-11/2-12/3-13/4-14/5-15/6-16/7-17/8-18/9-19/10-20/11-12/11-18/11-21/12-20/12-21/13-14/13-15/13-21/14-16/14-21/15-17/15-21/16-17/16-21/17-21/18-19/18-21/19-20/19-21/20-21",
        ),
        (
            "bipyridine",
            "C10H8N2/1-9/2-10/3-11/4-15/5-12/6-13/7-14/8-16/9-10/9-11/10-15/11-17/12-13/12-14/13-16/14-18/15-19/16-20/17-18/17-19/18-20",
        ),
        (
            "cf3alkyne",
            "C6H5F3O2/1-6/2-6/3-6/4-9/5-9/6-9/7-8/7-10/8-11/9-13/10-12/10-13/11-14/11-15/11-16",
        ),
    ],
)
def test_regression(m, expected_serialization):
    m = graph_from_file((Path(f"tests/molfiles/{m}/{m}.mol")))
    m_serialized = serialize_molecule(canonicalize_molecule(m, 0))
    assert m_serialized == expected_serialization
