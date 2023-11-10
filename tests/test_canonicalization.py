from tucan.canonicalization import (
    canonicalize_molecule,
    partition_molecule_by_attribute,
    refine_partitions,
)
from tucan.graph_attributes import ATOMIC_NUMBER
from tucan.serialization import serialize_molecule
from tucan.graph_utils import permute_molecule
from tucan.test_utils import permutation_invariance
from tucan.io import graph_from_file
from itertools import pairwise
import networkx as nx
import pytest


@pytest.mark.parametrize(
    "m, expected_partitions",
    [
        ("tnt", [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4]),
        ("favipiravir_1", [0, 1, 1, 2, 3, 4, 5, 5, 6, 7, 8, 8, 9, 10, 11]),
    ],
)
def test_partition_molecule_by_attribute(m, expected_partitions):
    m_partitioned = partition_molecule_by_attribute(
        graph_from_file(f"tests/molfiles/{m}/{m}.mol"), ATOMIC_NUMBER
    )
    partitions = sorted(nx.get_node_attributes(m_partitioned, "partition").values())

    assert partitions == expected_partitions


def test_partition_molecule_by_attribute_is_stable(m):
    m_partitioned = partition_molecule_by_attribute(m, ATOMIC_NUMBER)
    m_re_partitioned = partition_molecule_by_attribute(m_partitioned, ATOMIC_NUMBER)

    assert sorted(
        nx.get_node_attributes(m_partitioned, "partition").values()
    ) == sorted(nx.get_node_attributes(m_re_partitioned, "partition").values())


def test_refine_partitions(m):
    m_partitioned = partition_molecule_by_attribute(m, ATOMIC_NUMBER)
    m_refined = list(refine_partitions(m_partitioned))

    assert all(
        max_p_i < max_p_j
        for max_p_i, max_p_j in pairwise(
            (
                max(nx.get_node_attributes(m_rfnd, "partition").values())
                for m_rfnd in m_refined
            )
        )
    )


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
