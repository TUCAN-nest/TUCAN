
import networkx as nx

from chembl_webresource_client.new_client import new_client

from tucan.element_properties import ELEMENT_PROPS

def graph_from_chembl(id: int) -> nx.Graph:
    """Instantiate a NetworkX graph from a ChEMBL entry.
    Parameters
    ----------
    id: int
        Non-zero integer number pointing to a ChEMBL entry [1].
    Returns
    -------
    NetworkX Graph, error (is either True if correctly processed or False if problem was encountered)
    References
    ----------
    [1] https://github.com/chembl/chembl_webresource_client
    """

    if(id < 1):
        raise ChemblParserException(
            f'Invalid Compound ID (ID) "{cid}"'
        )

    molecule = new_client.molecule

    chembl_id = 'CHEMBL'+str(id)

    error = False
    try:
        structure = molecule.filter(chembl_id = chembl_id).only(['molecule_structures'])[0]['molecule_structures']
    except:
        error = True
        m = []

    if(error == False):

        try:
            molfile2000 = structure['molfile']
        except:
            error = True
            m = []

        if(error == False):
            lines = molfile2000.splitlines()

            atoms = int(lines[3][0:3])
            bonds = int(lines[3][4:6])

            m = nx.Graph()

            for a in range(0, atoms):
                line = lines[4+a]
                element = line[31:34]
                element = element.replace(" ", "")
                if(element == "D"):
                    element = "H"
                if(element == "T"):
                    element = "H"
                m.add_node(a+1)
                atom_props = { a+1: {
                                        "element_symbol": element,
                                        "atomic_number": ELEMENT_PROPS[element]["atomic_number"]
                                    }
                             }
                nx.set_node_attributes(m, atom_props)
            for b in range(0, bonds):
                line = lines[4+a+1+b]
                atom1 = int(line[0:3])
                atom2 = int(line[4:6])
                m.add_edge(atom1, atom2)

    return m, error

class ChemblParserException(Exception):
    pass
