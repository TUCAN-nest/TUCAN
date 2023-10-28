
import networkx as nx
import pubchempy as pcp

from tucan.element_properties import ELEMENT_PROPS

def graph_from_pubchem(cid: int) -> nx.Graph:
    """Instantiate a NetworkX graph from a PubChem entry.
    Parameters
    ----------
    cid: int
        Non-zero integer number pointing to a PubChem entry [1].
    Returns
    -------
    NetworkX Graph
    References
    ----------
    [1] https://pubchempy.readthedocs.io
    """

    if(cid < 1):
        raise PubChemParserException(
            f'Invalid Compound ID (CID) "{cid}"'
        )

    c = pcp.Compound.from_cid(cid)

    atoms = c.to_dict(properties=['atoms'])
    bonds = c.to_dict(properties=['bonds'])

    m = nx.Graph()

    for atom in atoms["atoms"]:
        keys = atom.keys()
        for k in keys:
            if(k == "aid"):
                atom1 = atom[k]
            elif(k == "number"):
                atomic_number = atom[k]
            elif(k == "element"):
                element_symbol = atom[k]
            elif(k == "x"):
                xcoord = atom[k]
            elif(k == "y"):
                ycoord = atom[k]
        zcoord = 0
        m.add_node(int(atom1))
        attrs = {atom1: {"node_label": atom1, "atomic_number": atomic_number, "partition": 0, "element_symbol": element_symbol, "element_color": (208,208,224),
                         "x_coord": float(xcoord), "y_coord": float(ycoord), "z_coord": float(zcoord)
                        }
                }
        nx.set_node_attributes(m, attrs)

    for bond in bonds["bonds"]:
        keys = bond.keys()
        for k in keys:
            if(k == "aid1"):
                atom1 = bond[k]
            elif(k == "aid2"):
                atom2 = bond[k]
            elif(k == "order"):
                bond_order = bond[k]
        m.add_edge(int(atom1), int(atom2))

    return m

class PubChemParserException(Exception):
    pass
