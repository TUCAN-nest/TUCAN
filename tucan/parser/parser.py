from antlr4 import InputStream, CommonTokenStream
from antlr4.error.ErrorListener import ErrorListener
from tucan.parser.tucanLexer import tucanLexer
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
