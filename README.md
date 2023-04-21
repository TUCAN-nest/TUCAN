![CI](https://github.com/TUCAN-nest/TUCAN/actions/workflows/ci.yml/badge.svg)


<img src="https://github.com/TUCAN-nest/TUCAN/blob/HEAD/logo.png" alt="logo" style="width:400px;"/>

TUCAN[^1] is a molecular identifier and descriptor for all domains of chemistry.

When you give TUCAN a molfile:

```
HOAc.mol
 OpenBabel02172217462D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.357235 0 0 0
M  V30 2 C 0.357235 0.4125 0 0
M  V30 3 H 1.07171 0 0 0
M  V30 4 H 0.357235 1.2375 0 0
M  V30 5 H 1.07171 0.825 0 0
M  V30 6 O -1.07171 0.4125 0 0
M  V30 7 O -0.357235 -0.825 0 0
M  V30 8 H 0.357235 -1.2375 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 1 6
M  V30 3 1 1 7
M  V30 4 1 2 3
M  V30 5 1 2 4
M  V30 6 1 2 5
M  V30 7 1 7 8
M  V30 END BOND
M  V30 END CTAB
M  END
```

it computes a serialized, canonical representation of your molecule:

```Python
from tucan.io.molfile_reader import graph_from_file
from tucan.canonicalization import canonicalize_molecule
from tucan.serialization import serialize_molecule

molecule = graph_from_file("path/to/your/molecule.mol")
canonical_molecule = canonicalize_molecule(molecule)
tucan_string = serialize_molecule(canonical_molecule)

print(tucan_string)
```

```bash
C2H4O2/(1-5)(2-5)(3-5)(4-8)(5-6)(6-7)(6-8)
```

Have a look at our [demo notebook](https://github.com/TUCAN-nest/TUCAN/blob/HEAD/docs/demo.ipynb) to see TUCAN at work.

You can explore TUCAN interactively on our [website](https://tucan-nest.github.io),

and find details on TUCAN's scientific background in our [paper](https://dx.doi.org/10.1186/s13321-022-00640-5).

# Installation
TUCAN requires a Python installation (>=3.10), and preferably you install the package in a virtual environment (e.g., [venv](https://docs.python.org/3.11/library/venv.html)).

You can install TUCAN locally using [pip](https://pip.pypa.io/en/stable/) by first [cloning the repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository), and then running
```
pip install .
```
in the root of the repository. Alternatively, you can skip cloning by running
```
pip install git+https://github.com/TUCAN-nest/TUCAN.git
```

In order to install an editable development version run
```
pip install -e .
```
The commands above will only install the minimal set of dependencies that is necessary to use TUCAN. If you want to install the additional development dependencies (testing, linting, support for notebooks), run

```
pip install .[dev,drawing]
```
or
```
pip install git+https://github.com/TUCAN-nest/TUCAN.git#egg=tucan[dev,drawing]
```

Have a look at the [pyproject.toml](https://github.com/TUCAN-nest/TUCAN/blob/HEAD/pyproject.toml) for details regarding the optional dependencies.

# Contributors
* Claudia Kellner (Universität Würzburg) - molecule test set
* Alexander Hoffmann (RWTH Aachen) - molecule test set
* Sonja Herres-Pawlis (RWTH Aachen) - concept
* Ulrich Schatzschneider (Universität Würzburg) - concept, documentation, and software development
* Jan C. Brammer (RWTH Aachen) - software development and documentation
* Frank Lange (IPB Halle) - software development

[^1]:TUCAN stands for **tu**ple **can**onicalization, since we're using 2-tuples to represent the bonds of our canonicalized molecules. Besides, [Wikipedia informs us](https://en.wikipedia.org/wiki/Toco_toucan) that in the tucán's native parlance, "_tu canonicalización_" expresses our hope that this software will help canonicalize _your_ molecules.
