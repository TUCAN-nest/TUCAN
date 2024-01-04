from tabulate import tabulate
import matplotlib.pyplot as plt
import networkx as nx
import plotly.graph_objects as go
import plotly.subplots as sp

from tucan.graph_attributes import (
    ATOMIC_NUMBER,
    INVARIANT_CODE,
    PARTITION,
)


def _draw_networkx_graph(m, highlight, labels, ax):
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


def _draw_networkx_graph_3d(m, highlight, labels, fig, col):
    coords = nx.kamada_kawai_layout(m, dim=3)
    highlight_colors = list(nx.get_node_attributes(m, highlight).values())
    labels = list(map(str, list(m.nodes))) if not labels else labels

    # Plotly requires separate node coordinates...
    x_nodes = [coords[key][0] for key in coords.keys()]
    y_nodes = [coords[key][1] for key in coords.keys()]
    z_nodes = [coords[key][2] for key in coords.keys()]
    # ...as well as tuples of node coordinated that define edges.
    x_edges = []
    y_edges = []
    z_edges = []
    for edge in m.edges():
        x_edges.extend([coords[edge[0]][0], coords[edge[1]][0], None])
        y_edges.extend([coords[edge[0]][1], coords[edge[1]][1], None])
        z_edges.extend([coords[edge[0]][2], coords[edge[1]][2], None])

    trace_edges = go.Scatter3d(
        x=x_edges,
        y=y_edges,
        z=z_edges,
        mode="lines",
        line=dict(color="black", width=2),
        hoverinfo="none",
    )

    trace_nodes = go.Scatter3d(
        x=x_nodes,
        y=y_nodes,
        z=z_nodes,
        mode="markers+text",
        marker=dict(
            symbol="circle",
            size=7,
            color=highlight_colors,
            colorscale="plotly3_r",  # turbo
            opacity=0.5,
        ),
        text=labels,
        hoverinfo="none",
        textfont=dict(size=7.5),
    )

    fig.add_traces([trace_edges, trace_nodes], rows=1, cols=col)


def draw_molecules(
    m_list, caption_list, labels=None, highlight=ATOMIC_NUMBER, title="", dim=2
):
    """Draw molecule(s).

    Parameters
    ----------
    highlight: str
        Color atoms by "atomic_number" (default) or "partition".
    dim: int
        Plot in "2" (default) or "3" dimensions.
    """
    if highlight not in [ATOMIC_NUMBER, PARTITION]:
        print(
            f"Please select one of {{'{PARTITION}', '{ATOMIC_NUMBER}'}} for `highlight`."
        )
        return
    n_molecules = len(m_list)

    if dim == 2:
        fig = plt.figure(figsize=(n_molecules * 6, 6))
        fig.suptitle(title)
        for i, m in enumerate(m_list):
            ax = fig.add_subplot(1, n_molecules, i + 1, title=caption_list[i])
            _draw_networkx_graph(m, highlight, labels, ax)
    elif dim == 3:
        fig = sp.make_subplots(
            rows=1,
            cols=n_molecules,
            subplot_titles=caption_list,
            specs=[[{"is_3d": True}] * n_molecules],
            horizontal_spacing=0.1 / n_molecules,
        )
        for i, m in enumerate(m_list):
            _draw_networkx_graph_3d(m, highlight, labels, fig, i + 1)
        fig.update_layout(showlegend=False, title=title)
        fig.update_scenes(
            patch=dict(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False),
                camera_eye=dict(x=1, y=1, z=1),
            ),
        )
        fig.show()


def print_molecule(m, caption=""):
    print(caption)
    table = []
    for atom in sorted(list(m.nodes)):
        invariant_code = m.nodes[atom][INVARIANT_CODE]
        partition = m.nodes[atom][PARTITION]
        neighbors = [
            (n, m.nodes[n][INVARIANT_CODE], m.nodes[n][PARTITION])
            for n in m.neighbors(atom)
        ]
        neighbors = sorted(neighbors, key=lambda x: x[2], reverse=True)
        neighbors = ", ".join([f"({i}, {c}, {p})" for i, c, p in neighbors])
        table.append([atom, invariant_code, partition, neighbors])
    print(
        tabulate(
            table,
            tablefmt="simple",
            colalign=["left"] * 4,
            headers=[
                "label",
                "invariants",
                "partition",
                "neighbors (label, invariants, partition)",
            ],
        )
    )
