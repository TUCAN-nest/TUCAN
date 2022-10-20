from tucan.io import graph_from_file
from tucan.serialization import serialize_molecule
from tucan.canonicalization import canonicalize_molecule
from pathlib import Path
import pytest


@pytest.mark.parametrize(
    "m, expected_serialization",
    [
        (
            "Bicyclo[5.1.0]octa-1(7)-ene-8-one",
            "C8H10O/(1-11)(2-11)(3-12)(4-12)(5-13)(6-13)(7-14)(8-14)(9-15)(10-15)(11-12)(11-13)(12-14)(13-15)(14-16)(15-17)(16-17)(16-18)(17-18)(18-19)",
        ),
        (
            "ferrocene",
            "C10H10Fe/(1-11)(2-12)(3-13)(4-14)(5-15)(6-16)(7-17)(8-18)(9-19)(10-20)(11-12)(11-13)(11-21)(12-14)(12-21)(13-15)(13-21)(14-15)(14-21)(15-21)(16-17)(16-18)(16-21)(17-19)(17-21)(18-20)(18-21)(19-20)(19-21)(20-21)",
        ),
        (
            "bipyridine",
            "C10H8N2/(1-9)(2-10)(3-11)(4-12)(5-13)(6-14)(7-15)(8-16)(9-11)(9-13)(10-12)(10-14)(11-15)(12-16)(13-17)(14-18)(15-19)(16-20)(17-18)(17-19)(18-20)",
        ),
        (
            "cf3alkyne",
            "C6H5F3O2/(1-6)(2-6)(3-6)(4-9)(5-9)(6-9)(7-8)(7-10)(8-11)(9-13)(10-12)(10-13)(11-14)(11-15)(11-16)",
        ),
        (
            "tnt",
            "C6H3N3O6/(1-4)(2-5)(3-6)(4-7)(4-8)(5-7)(5-9)(6-8)(6-9)(7-10)(8-11)(9-12)(10-13)(10-14)(11-15)(11-16)(12-17)(12-18)",
        ),
        (
            "TEMPO",
            "C9H18NO/(1-19)(2-19)(3-19)(4-20)(5-20)(6-20)(7-21)(8-21)(9-21)(10-22)(11-22)(12-22)(13-23)(14-23)(15-24)(16-24)(17-25)(18-25)(19-26)(20-26)(21-27)(22-27)(23-24)(23-25)(24-26)(25-27)(26-28)(27-28)(28-29)/(29:rad=2)",
        ),
        (
            "water-d1_1",
            "H2O/(1-3)(2-3)/(2:mass=2)",
        ),
        (
            "water-d1_3",
            "H2O/(1-3)(2-3)/(2:mass=2)",
        ),
        (
            "water-d2",
            "H2O/(1-3)(2-3)/(1:mass=2)(2:mass=2)",
        ),
        (
            "water-t2",
            "H2O/(1-3)(2-3)/(1:mass=3)(2:mass=3)",
        ),
    ],
)
def test_regression(m, expected_serialization):
    m = graph_from_file((Path(f"tests/molfiles/{m}/{m}.mol")))
    m_serialized = serialize_molecule(canonicalize_molecule(m))
    assert m_serialized == expected_serialization
