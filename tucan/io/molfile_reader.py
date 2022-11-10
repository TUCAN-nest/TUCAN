import networkx as nx
from pathlib import Path
from tucan.io.exception import MolfileParserException
from tucan.io.molfile_v2000_reader import graph_from_molfile_v2000
from tucan.io.molfile_v3000_reader import graph_from_molfile_v3000


def graph_from_file(filepath: str) -> nx.Graph:
    """Instantiate a NetworkX graph from an MDL molfile. The parser supports both
    V3000 and V2000 connection tables.

    Parameters
    ----------
    filepath: str
        Path pointing to a file containing an MDL (now BIOVIA) connection
        table [1,2].

    Returns
    -------
    NetworkX Graph

    References
    ----------
    [1] https://en.wikipedia.org/wiki/Chemical_table_file
    [2] https://discover.3ds.com/sites/default/files/2020-08/biovia_ctfileformats_2020.pdf
    """
    filepath = Path(filepath)
    if filepath.suffix != ".mol":
        raise IOError(f"The file must be in '.mol' format, not {filepath.suffix}.")
    with open(filepath) as file:
        filecontent = file.read()

    return graph_from_molfile_text(filecontent)


def graph_from_molfile_text(molfile: str) -> nx.Graph:
    lines = molfile.splitlines()

    molfile_version = lines[3].rstrip().split(" ")[-1]
    if molfile_version == "V3000":
        return graph_from_molfile_v3000(lines)
    elif molfile_version == "V2000":
        return graph_from_molfile_v2000(lines)
    else:
        raise MolfileParserException(
            f'Unsupported Molfile version "{molfile_version}"'
        )
