from antlr4 import InputStream, CommonTokenStream
from antlr4.error.ErrorListener import ErrorListener
from antlr4.tree.Tree import ParseTreeWalker
from tucan.element_properties import ELEMENT_PROPS
from tucan.parser.tucanLexer import tucanLexer
from tucan.parser.tucanListener import tucanListener
from tucan.parser.tucanParser import tucanParser


def parse_tucan(s):
    parser = _prepare_parser(s)
    tree = parser.tucan()
    # TODO: extract graph


def _prepare_parser(s):
    stream = InputStream(s)
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

        for _ in range(count):
            self._atoms.append(atom_props)

    def _add_bond(self, index1, index2):
        self._bonds.append((index1 - 1, index2 - 1))

    def _add_node_property(self, node_index, key, value):
        props_for_node = self._node_attributes.setdefault(node_index, {})

        if key in props_for_node:
            raise TucanParserException(
                f'Atom {node_index}: Property "{key}" was already defined.'
            )
        props_for_node[key] = value


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
