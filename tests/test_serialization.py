from tucan.serialization import serialize_molecule
from tucan.canonicalization import canonicalize_molecule


def test_regression(m, snapshot):
    """
    Generate TUCAN strings for all test molfiles and match them against the
    snapshot file tests/__snapshots__/test_serialization.ambr.

    Run `pytest --snapshot-update` to update the snapshot file.
    """
    m_serialized = serialize_molecule(canonicalize_molecule(m))
    assert m_serialized == snapshot
