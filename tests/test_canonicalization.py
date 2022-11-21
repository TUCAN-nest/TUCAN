from tucan.canonicalization import canonicalize_molecule
from tucan.serialization import serialize_molecule
from tucan.graph_utils import permute_molecule
from tucan.tests import permutation_invariance
import networkx as nx
import pytest


def test_permutation(m):
    # Enforce permutation for graphs with at least 2 edges that aren't fully connected (i.e., complete).
    if m.number_of_edges() <= 1:
        pytest.skip("Skipping graph with less than two edges.")
    if nx.density(m) == 1:
        pytest.skip("Skipping fully connected graph.")

    identity = {label: label for label in m}
    nx.set_node_attributes(m, identity, "original_label")

    m_permu = permute_molecule(m, random_seed=0.42)
    original_labels = [data for _, data in m.nodes("original_label")]
    permuted_labels = [data for _, data in m_permu.nodes("original_label")]

    # checks that iteration order of the graph nodes has indeed changed
    assert original_labels != permuted_labels

    assert m.edges != m_permu.edges


def test_permutation_invariance(m):
    permutation_invariance(m)


@pytest.mark.parametrize(
    "excludes",
    [
        [
            "water-d1_3",  # duplicate of water-d1_1
        ]
    ],
)
def test_bijection(excludes):
    """Eineindeutigkeit."""
    serializations = set()
    for f in pytest.testset:
        if f.stem in excludes:
            continue

        m = pytest.filereader(f)
        m_serialized = serialize_molecule(canonicalize_molecule(m))
        assert m_serialized not in serializations, f"duplicate: {f.stem}"
        serializations.add(m_serialized)
