from tabulate import tabulate
import matplotlib.pyplot as plt
import networkx as nx


def _rdkit_to_nx(m):
    """Convert an RDKit molecule to a NetworkX graph with node properties
    "atomic_number" and "partition"."""
    nx_graph = nx.Graph()
    nx_graph.add_nodes_from([a.GetIdx() for a in m.GetAtoms()])
    for b in m.GetBonds():
        nx_graph.add_edge(b.GetBeginAtom().GetIdx(), b.GetEndAtom().GetIdx())
    atomic_numbers = [a.GetAtomicNum() for a in m.GetAtoms()]
    nx.set_node_attributes(nx_graph, {i:p for i, p in enumerate(atomic_numbers)}, "atomic_number")
    partitions = [a.GetIntProp("partition") if a.HasProp("partition") else 0 for a in m.GetAtoms()]
    nx.set_node_attributes(nx_graph, {i:p for i, p in enumerate(partitions)}, "partition")
    return nx_graph

def _draw_networkx_graph(m, ax, highlight):
  nx_graph = _rdkit_to_nx(m) if type(m) != nx.Graph else m
  highlight_colors = list(nx.get_node_attributes(nx_graph, highlight).values())
  node_size = 1 / nx_graph.order() * 10000
  nx.draw_kamada_kawai(nx_graph, node_color=highlight_colors, node_size=node_size,
                       cmap="rainbow", alpha=.5, with_labels=True, font_weight="heavy", ax=ax)

def draw_molecules(m_list, caption_list, highlight="atomic_number", title=""):
  """`highlight`: color atoms by "atomic_number" (default) or "partition"."""
  if highlight not in ["atomic_number", "partition"]:
    print("Please select one of {'partition', 'atomic_number'} for `highlight`.")
    return
  n_molecules = len(m_list)
  fig = plt.figure(figsize=(n_molecules * 6, 6))
  fig.suptitle(title)
  for i, m in enumerate(m_list):
    ax = fig.add_subplot(1, n_molecules, i + 1, title=caption_list[i])
    _draw_networkx_graph(m, ax, highlight)

def print_molecule(m, caption=""):
    nx_graph = _rdkit_to_nx(m) if type(m) != nx.Graph else m
    print(caption)
    table = []
    for atom in nx_graph.nodes():
        fingerprint = nx_graph.nodes[atom]["fingerprint"]
        partition = nx_graph.nodes[atom]["partition"]
        neighbors = [(n, nx_graph.nodes[n]["fingerprint"], nx_graph.nodes[n]["partition"])
                     for n in nx_graph.neighbors(atom)]
        table.append([atom, fingerprint, partition, neighbors])
    print(tabulate(table, tablefmt="fancy_grid",
                   headers=["index", "fingerprint", "partition",
                            "neighbors (index, fingerprint, partition)"]))
