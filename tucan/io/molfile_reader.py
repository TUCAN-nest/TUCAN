import networkx as nx
from pathlib import Path
from tucan.io.exception import MolfileParserException
from tucan.io.molfile_v2000_reader import graph_props_from_molfile_v2000
from tucan.io.molfile_v3000_reader import graph_props_from_molfile_v3000


def graph_from_file(filepath: str) -> nx.Graph:
    """Instantiate a NetworkX graph from an MDL molfile as specified by [1]. The
    parser supports both V3000 and V2000 connection tables (CTABs).

    Parameters
    ----------
    filepath: str
        Path pointing to a molfile (.mol file extension)

    Returns
    -------
    NetworkX Graph

    References
    ----------
    [1] https://discover.3ds.com/sites/default/files/2020-08/biovia_ctfileformats_2020.pdf
    """
    filepath_object = Path(filepath)
    if filepath_object.suffix != ".mol":
        raise IOError(f"The file must be in '.mol' format, not {filepath_object.suffix}.")
    with open(filepath_object) as file:
        filecontent = file.read()

    return graph_from_molfile_text(filecontent)


def graph_from_molfile_text(molfile: str) -> nx.Graph:
    """Instantiate a NetworkX graph from an MDL molfile as specified by [1]. The
    parser supports both V3000 and V2000 connection tables (CTABs).

    Parameters
    ----------
    molfile: str
        the molfile as string

    Returns
    -------
    NetworkX Graph

    References
    ----------
    [1] https://discover.3ds.com/sites/default/files/2020-08/biovia_ctfileformats_2020.pdf
    """
    lines = molfile.splitlines()

    molfile_version = lines[3].rstrip().split(" ")[-1]
    if molfile_version == "V3000":
        atom_props, bond_props = graph_props_from_molfile_v3000(lines)
    elif molfile_version == "V2000":
        atom_props, bond_props = graph_props_from_molfile_v2000(lines)
    else:
        raise MolfileParserException(f'Unsupported Molfile version "{molfile_version}"')

    graph = nx.Graph()
    graph.add_nodes_from(list(atom_props.keys()))
    nx.set_node_attributes(graph, atom_props)
    graph.add_edges_from(list(bond_props.keys()))
    nx.set_edge_attributes(graph, bond_props)

    return nx.convert_node_labels_to_integers(graph)
