# Parsing TUCAN strings

[parser.py](parser.py) contains a recursive-descent parser.
The TUCAN grammar is documented in [tucan.g4](tucan.g4).

When changing the TUCAN grammar, manually update `tucan.g4` and (if necessary) the corresponding tests under `tests/parser`,
and then instruct an agent to update `parser.py` according to those constraints.
