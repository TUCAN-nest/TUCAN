from tabulate import tabulate
import matplotlib.pyplot as plt
import networkx as nx


def _draw_networkx_graph(m, ax, labels, highlight):
    highlight_colors = list(nx.get_node_attributes(m, highlight).values())
    node_size = 1 / m.order() * 10000
    nx.draw_kamada_kawai(
        m,
        node_color=highlight_colors,
        node_size=node_size,
        cmap="rainbow",
        alpha=0.5,
        with_labels=True,
        labels=labels,
        font_weight="heavy",
        ax=ax,
    )


def draw_molecules(
    m_list, caption_list, labels=None, highlight="atomic_number", title=""
):
    """Draw molecule(s).

    Parameters
    ----------
    highlight: str
        Color atoms by "atomic_number" (default) or "partition".
    """
    if highlight not in ["atomic_number", "partition"]:
        print("Please select one of {'partition', 'atomic_number'} for `highlight`.")
        return
    n_molecules = len(m_list)
    fig = plt.figure(figsize=(n_molecules * 6, 6))
    fig.suptitle(title)
    for i, m in enumerate(m_list):
        ax = fig.add_subplot(1, n_molecules, i + 1, title=caption_list[i])
        _draw_networkx_graph(m, ax, labels, highlight)


def print_molecule(m, caption=""):
    print(caption)
    table = []
    for atom in sorted(m.nodes()):
        invariant_code = m.nodes[atom]["invariant_code"]
        partition = m.nodes[atom]["partition"]
        neighbors = [
            (n, m.nodes[n]["invariant_code"], m.nodes[n]["partition"])
            for n in m.neighbors(atom)
        ]
        neighbors = sorted(neighbors, key=lambda x: x[1], reverse=True)
        neighbors = ", ".join([f"({i}, {c}, {p})" for i, c, p in neighbors])
        table.append([atom, invariant_code, partition, neighbors])
    print(
        tabulate(
            table,
            tablefmt="simple",
            colalign=["left"] * 4,
            headers=[
                "label",
                "invariant-code",
                "partition",
                "neighbors (label, invariant-code, partition)",
            ],
        )
    )
