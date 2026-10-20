from __future__ import annotations

import networkx as nx
from typing import Any, Final
from tucan.element_attributes import ELEMENT_ATTRS, element_symbols
from tucan.graph_attributes import (
    ATOMIC_NUMBER,
    ELEMENT_SYMBOL,
    PARTITION,
)
from tucan.graph_utils import graph_from_molecule
from tucan.serialization import _SERIALIZER_NODE_ATTRIBUTE_MAPPING


def graph_from_tucan(tucan: str) -> nx.Graph:
    """Instantiate a NetworkX graph from a TUCAN string.

    Parameters
    ----------
    tucan: str
        TUCAN string to be deserialized.

    Returns
    -------
    NetworkX Graph
    """
    return _Parser(tucan).parse_tucan().to_graph()


_DESERIALIZER_NODE_ATTRIBUTE_MAPPING: Final[dict[str, str]] = {
    v: k for k, v in _SERIALIZER_NODE_ATTRIBUTE_MAPPING.items()
}


class _GraphBuilder:
    def __init__(self):
        self._atoms = []
        self._bonds = []
        self._node_attributes = {}  # becomes a dictionary of dictionaries

    def add_sum_formula(self, symbols_and_counts):
        for symbol, count in symbols_and_counts:
            self._add_atoms(symbol, count)

    def add_tuple(self, index1, index2):
        if index1 == index2:
            raise TucanParserException(
                f'Error in tuple "({index1}-{index2})": Self-loops are not allowed.'
            )
        self._add_bond(index1, index2)

    def add_node_property(self, node_index, key, value):
        self._add_node_attribute(node_index, key, value)

    def _add_atoms(self, element, count):
        atom_attrs = {
            ELEMENT_SYMBOL: element,
            ATOMIC_NUMBER: ELEMENT_ATTRS[element][ATOMIC_NUMBER],
            PARTITION: 0,
        }

        self._atoms.extend([atom_attrs.copy() for _ in range(count)])

    def _add_bond(self, index1, index2):
        self._bonds.append((index1 - 1, index2 - 1))

    def _add_node_attribute(self, node_index, key, value):
        attrs_for_node = self._node_attributes.setdefault(node_index - 1, {})
        attr_key = _DESERIALIZER_NODE_ATTRIBUTE_MAPPING[key]

        if attr_key in attrs_for_node:
            raise TucanParserException(
                f'Atom {node_index}: Attribute "{key}" was already defined.'
            )
        attrs_for_node[attr_key] = value

    def to_graph(self) -> nx.Graph:
        # node index validation
        for i1, i2 in self._bonds:
            self._validate_atom_index(i1)
            self._validate_atom_index(i2)

        sorted_atoms = sorted(self._atoms, key=lambda a: a[ATOMIC_NUMBER])

        # dict of dict (atom_index -> dict of atom attributes)
        atoms_dict: dict[int, dict[str, Any]] = {
            i: sorted_atoms[i] for i in range(len(sorted_atoms))
        }

        # join in additional atom attributes
        for index, attrs in self._node_attributes.items():
            self._validate_atom_index(index)

            atom_attrs = atoms_dict[index]
            atom_attrs.update(attrs)

        bonds_dict: dict[tuple[int, int], dict] = {bond: {} for bond in self._bonds}

        return graph_from_molecule(atoms_dict, bonds_dict)

    def _validate_atom_index(self, index):
        if index >= len(self._atoms):
            raise TucanParserException(f"Atom with index {index + 1} does not exist.")


_ELEMENT_SYMBOLS: Final[frozenset[str]] = frozenset(element_symbols)
_NODE_PROPERTY_KEYS: Final[tuple[str, ...]] = ("mass", "rad")
_VERSION_PREFIX: Final[str] = "TUCANv"


class _Parser:
    """Recursive-descent parser for TUCAN strings (see tucan.g4 for the grammar)."""

    def __init__(self, text: str):
        self._text = text
        self._pos = 0
        self._builder = _GraphBuilder()

    # entry points
    def parse_tucan(self) -> _GraphBuilder:
        self._version()
        self._expect("/")
        self._sum_formula()
        self._expect("/")
        self._tuples()
        if self._accept("/"):
            self._node_attributes()
        self._expect_end()
        return self._builder

    # error handling
    def _error(self, msg: str, pos: int | None = None, length: int = 1):
        pos = self._pos if pos is None else pos
        line_no = self._text.count("\n", 0, pos) + 1
        line_start = self._text.rfind("\n", 0, pos) + 1
        line_end = self._text.find("\n", pos)
        if line_end == -1:
            line_end = len(self._text)
        column = pos - line_start
        marker = f"{column * ' '}{max(length, 1) * '^'}"
        raise TucanParserException(
            f"line {line_no}:{column} {msg}\n{self._text[line_start:line_end]}\n{marker}"
        )

    def _found(self) -> str:
        return repr(self._text[self._pos]) if self._pos < len(self._text) else "<EOF>"

    # token helpers
    def _peek(self) -> str:
        return self._text[self._pos : self._pos + 1]

    def _accept(self, literal: str) -> bool:
        if self._text.startswith(literal, self._pos):
            self._pos += len(literal)
            return True
        return False

    def _expect(self, literal: str) -> None:
        if not self._accept(literal):
            self._error(f"expected '{literal}' but found {self._found()}")

    def _expect_end(self) -> None:
        if self._pos < len(self._text):
            self._error(f"expected end of input but found {self._found()}")

    def _number(self, minimum: int) -> int:
        start = self._pos
        if self._peek() == "0":
            self._pos += 1
        elif self._peek() and self._peek() in "123456789":
            while self._peek().isdecimal() and self._peek().isascii():
                self._pos += 1
        else:
            self._error(f"expected a number but found {self._found()}")
        value = int(self._text[start : self._pos])
        if value < minimum:
            self._error(
                f"expected a number >= {minimum} but found '{value}'",
                start,
                self._pos - start,
            )
        return value

    # grammar rules
    def _version(self) -> None:
        self._expect(_VERSION_PREFIX)
        self._number(0)
        self._expect(".")
        self._number(0)
        self._expect(".")
        self._number(0)

    def _sum_formula(self) -> None:
        entries = []
        while self._peek().isascii() and self._peek().isupper():
            start = self._pos
            symbol = self._element_symbol()
            count = self._number(2) if self._peek().isdecimal() else 1
            entries.append((start, symbol, count))
        self._validate_hill_order(entries)
        self._builder.add_sum_formula((symbol, count) for _, symbol, count in entries)

    def _element_symbol(self) -> str:
        for length in (2, 1):
            symbol = self._text[self._pos : self._pos + length]
            if len(symbol) == length and symbol in _ELEMENT_SYMBOLS:
                self._pos += length
                return symbol
        self._error(f"unknown element symbol at {self._found()}")
        raise AssertionError  # unreachable

    def _validate_hill_order(self, entries) -> None:
        remaining = list(entries)
        if remaining and remaining[0][1] == "C":
            remaining.pop(0)
            if remaining and remaining[0][1] == "H":
                remaining.pop(0)
            forbidden = {"C", "H"}
        else:
            forbidden = {"C"}

        previous = ""
        for start, symbol, _ in remaining:
            if symbol in forbidden or symbol <= previous:
                self._error(
                    f"element '{symbol}' violates the Hill order", start, len(symbol)
                )
            previous = symbol

    def _tuples(self) -> None:
        while self._accept("("):
            index1 = self._number(1)
            self._expect("-")
            index2 = self._number(1)
            self._expect(")")
            self._builder.add_tuple(index1, index2)

    def _node_attributes(self) -> None:
        while self._accept("("):
            node_index = self._number(1)
            self._expect(":")
            self._node_property(node_index)
            while self._accept(","):
                self._node_property(node_index)
            self._expect(")")

    def _node_property(self, node_index: int) -> None:
        key = next((k for k in _NODE_PROPERTY_KEYS if self._accept(k)), None)
        if key is None:
            self._error(f"expected 'mass' or 'rad' but found {self._found()}")
        self._expect("=")
        value = self._number(1)
        self._builder.add_node_property(node_index, key, value)


class TucanParserException(Exception):
    pass
