from tucan.serialization import serialize_molecule
from tucan.canonicalization import canonicalize_molecule


def test_regression(m, snapshot):
    m_serialized = serialize_molecule(canonicalize_molecule(m))
    snapshot.assert_match(m_serialized)
