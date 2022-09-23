import networkx as nx
from antlr4 import InputStream, CommonTokenStream
from antlr4.error.ErrorListener import ErrorListener
from antlr4.tree.Tree import ParseTreeWalker
from tucan.element_properties import ELEMENT_PROPS
from tucan.parser.tucanLexer import tucanLexer
from tucan.parser.tucanListener import tucanListener
from tucan.parser.tucanParser import tucanParser


def parse_tucan(tucan: str) -> nx.Graph:
    parser = _prepare_parser(tucan)
    tree = parser.tucan()
    listener = _walk_tree(tree)
    return listener.to_graph()


def _prepare_parser(to_parse: str) -> nx.Graph:
    stream = InputStream(to_parse)
    lexer = tucanLexer(stream)
    token_stream = CommonTokenStream(lexer)
    parser = tucanParser(token_stream)

    # register our own error listeners
    lexer.removeErrorListeners()
    lexer.addErrorListener(LexerErrorListener())
    parser.removeErrorListeners()
    parser.addErrorListener(ParserErrorListener())

    return parser


def _walk_tree(tree):
    walker = ParseTreeWalker()
    listener = TucanListenerImpl()
    walker.walk(listener, tree)
    return listener


class TucanListenerImpl(tucanListener):
    def __init__(self):
        self._atoms = []
        self._bonds = []
        self._node_attributes = {}  # becomes a dictionary of dictionaries

    def enterWith_carbon(self, ctx: tucanParser.With_carbonContext):
        self._parse_sum_formula(ctx)

    def enterWithout_carbon(self, ctx: tucanParser.Without_carbonContext):
        self._parse_sum_formula(ctx)

    def enterTuple(self, ctx: tucanParser.TupleContext):
        index1 = int(ctx.node_index(0).getText())
        index2 = int(ctx.node_index(1).getText())
        if index1 == index2:
            raise TucanParserException(
                f'Error in tuple "{ctx.getText()}": Self-loops are not allowed.'
            )
        self._add_bond(index1, index2)

    def enterNode_property(self, ctx: tucanParser.Node_propertyContext):
        node_index = int(ctx.parentCtx.node_index().getText())
        key = ctx.node_property_key().getText()
        value = int(ctx.node_property_value().getText())
        self._add_node_property(node_index, key, value)

    def _parse_sum_formula(self, formula_ctx):
        if formula_ctx.getChildCount() == 0:
            return

        for symbol_count_tuple in formula_ctx.children:
            symbol = symbol_count_tuple.getChild(0).getText()
            count = 1
            if symbol_count_tuple.getChildCount() > 1:
                count = int(symbol_count_tuple.getChild(1).getText())
            self._add_atoms(symbol, count)

    def _add_atoms(self, element, count):
        atom_props = {
            "element_symbol": element,
            "atomic_number": ELEMENT_PROPS[element]["atomic_number"],
            "partition": 0,
        }

        self._atoms.extend([atom_props.copy() for _ in range(count)])

    def _add_bond(self, index1, index2):
        self._bonds.append((index1 - 1, index2 - 1))

    def _add_node_property(self, node_index, key, value):
        props_for_node = self._node_attributes.setdefault(node_index - 1, {})

        if key in props_for_node:
            raise TucanParserException(
                f'Atom {node_index}: Property "{key}" was already defined.'
            )
        props_for_node[key] = value

    def to_graph(self) -> nx.Graph:
        # node index validation
        for i1, i2 in self._bonds:
            self._validate_atom_index(i1)
            self._validate_atom_index(i2)

        sorted_atoms = sorted(self._atoms, key=lambda a: a["atomic_number"])

        # dict of dict (atom_index -> dict of atom properties)
        atoms_dict = {i: sorted_atoms[i] for i in range(len(sorted_atoms))}

        # join in additional atom properties
        for index, props in self._node_attributes.items():
            self._validate_atom_index(index)

            atom_props = atoms_dict[index]
            atom_props.update(props)

        # construct graph
        graph = nx.Graph()
        graph.add_nodes_from(list(atoms_dict.keys()))
        nx.set_node_attributes(graph, atoms_dict)
        graph.add_edges_from(self._bonds)

        return graph

    def _validate_atom_index(self, index):
        if index >= len(self._atoms):
            raise TucanParserException(f"Atom with index {index + 1} does not exist.")


class RaisingErrorListener(ErrorListener):
    def syntaxError(self, recognizer, offending_symbol, line, column, msg, e):
        marked_error_location = self._underline_error(
            recognizer, offending_symbol, line, column
        )
        error_str = f"line {line}:{column} {msg}\n{marked_error_location}"
        raise TucanParserException(error_str)

    def _underline_error(self, recognizer, offending_symbol, line, column):
        pass


# The algorithm of _underline_error was adopted from Terence Parr's book "The Definitive ANTLR 4 Reference", page 156.
class LexerErrorListener(RaisingErrorListener):
    def _underline_error(self, recognizer, offending_symbol, line, column):
        lexer_input = str(recognizer.inputStream)
        error_line = lexer_input.split("\n")[line - 1]

        return f"{error_line}\n{column * ' '}^"


class ParserErrorListener(RaisingErrorListener):
    def _underline_error(self, recognizer, offending_symbol, line, column):
        tokens = recognizer.getInputStream()
        parser_input = str(tokens.tokenSource.inputStream)
        error_line = parser_input.split("\n")[line - 1]
        start = offending_symbol.start
        stop = offending_symbol.stop

        return f"{error_line}\n{column * ' '}{(stop - start + 1) * '^'}"


class TucanParserException(Exception):
    pass
