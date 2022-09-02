![CI](https://github.com/TUCAN-nest/TUCAN/actions/workflows/ci.yml/badge.svg)


<img src="https://github.com/TUCAN-nest/TUCAN/blob/HEAD/logo.png" alt="logo" style="width:400px;"/>

TUCAN[^1] - A molecular identifier and descriptor for all domains of chemistry

Have a look at our [demo notebook](https://github.com/TUCAN-nest/TUCAN/blob/HEAD/docs/demo.ipynb) to see TUCAN at work.

# Installation
TUCAN requires a Python installation (>=3.9), and preferably you install the package in a virtual environment (e.g., [venv](https://docs.python.org/3.10/library/venv.html)).

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

[^1]:TUCAN stands for **tu**ple **can**onicalization, since we're using 2-tuples to represent the bonds of our canonicalized molecules. Besides, [Wikipedia informs us](https://en.wikipedia.org/wiki/Toco_toucan) that in the tucán's native parlance, "_tu canonicalización_" expresses our hope that this software will help canonicalize _your_ molecules.
