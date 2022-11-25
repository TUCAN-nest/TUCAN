import networkx as nx
import random
from tucan.canonicalization import canonicalize_molecule
from tucan.graph_utils import permute_molecule
from tucan.parser.parser import graph_from_tucan
from tucan.serialization import serialize_molecule


def permutation_invariance(m: nx.Graph, n_runs=10, random_seed=random.random()):
    """Tests that different permutations of the same molecular graph yield one
    and the same TUCAN string."""
    m_canon = canonicalize_molecule(m)
    m_serialized = serialize_molecule(m_canon)
    random.seed(random_seed)
    for _ in range(n_runs):
        permutation_seed = random.random()
        m_permu = permute_molecule(m, random_seed=permutation_seed)
        m_permu_canon = canonicalize_molecule(m_permu)
        m_permu_serialized = serialize_molecule(m_permu_canon)
        assert m_serialized == m_permu_serialized


def roundtrip_graph_tucan_graph_tucan_graph(m: nx.Graph):
    m_serialized = serialize_molecule(canonicalize_molecule(m))
    g1 = graph_from_tucan(m_serialized)
    g1_serialized = serialize_molecule(canonicalize_molecule(g1.copy()))
    g2 = graph_from_tucan(g1_serialized)

    assert m_serialized == g1_serialized
    assert dict(g1.nodes(data=True)) == dict(g2.nodes(data=True))
    assert list(g1.edges) == list(g2.edges)
