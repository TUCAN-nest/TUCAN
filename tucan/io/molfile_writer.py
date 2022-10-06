from datetime import datetime
import networkx as nx
import tucan

_prog_name = f'TUCAN{tucan.__version__.replace(".", "")[0:3]: <3}'


def graph_to_molfile(graph: nx.Graph, calc_coordinates=False) -> str:
    """Generate an MDL V3000 Molfile from the given graph.

    Parameters
    ----------
    graph
        NetworkX Graph
    calc_coordinates:
        (optional) (re-)calculate atom positions

    Returns
    -------
    MDL V3000 Molfile
    """
    lines = []

    _add_header(lines)
    _add_v30_line(lines, "BEGIN CTAB")
    _add_v30_line(
        lines, f"COUNTS {graph.number_of_nodes()} {graph.number_of_edges()} 0 0 0"
    )
    _add_atom_block(lines, graph, calc_coordinates)
    _add_bond_block(lines, graph)
    _add_v30_line(lines, "END CTAB")
    lines.append(f"M  END")

    return "\n".join(lines)


def _add_header(lines: list[str]):
    # molecule name
    lines.append("")

    # IIPPPPPPPPMMDDYYHHmmddSSssssssssssEEEEEEEEEEEERRRRRR
    # A2<--A8--><---A10-->A2I2<--F10.5-><---F12.5--><-I6->
    # Fill until here:     ^
    lines.append(f"  {_prog_name}{datetime.now().strftime('%m%d%y%H%M')}3D")

    # comments line
    lines.append("")

    # CTAB version
    lines.append("  0  0  0     0  0            999 V3000")


def _add_v30_line(lines: list[str], line: str):
    # The length limit of a line is 80 characters. We include '\n' in this count.
    # "M  V30 " and '\n' take 8 chars, plus one char for '-' if line wrapping occurs.
    while True:
        if len(line) <= 72:
            lines.append(f"M  V30 {line}")
            break

        left, line = line[:71], line[71:]
        lines.append(f"M  V30 {left}-")


def _add_atom_block(lines: list[str], graph: nx.Graph, calc_coordinates: bool):
    if calc_coordinates:
        coords = nx.kamada_kawai_layout(graph, dim=2)

    _add_v30_line(lines, "BEGIN ATOM")

    for index, attrs in graph.nodes(data=True):
        x = coords[index][0] if calc_coordinates else attrs.get("x_coord", 0)
        y = coords[index][1] if calc_coordinates else attrs.get("y_coord", 0)
        z = 0 if calc_coordinates else attrs.get("z_coord", 0)

        charge = f" CHG={chg}" if (chg := attrs.get("chg")) and -15 <= chg <= 15 else ""
        radical = f" RAD={rad}" if (rad := attrs.get("rad")) and 0 < rad <= 3 else ""
        atomic_mass = (
            f" MASS={mass}" if (mass := attrs.get("mass")) and mass > 0 else ""
        )

        _add_v30_line(
            lines,
            f"{index + 1} {attrs['element_symbol']} {x:.6g} {y:.6g} {z:.6g} 0{charge}{radical}{atomic_mass}",
        )

    _add_v30_line(lines, "END ATOM")


def _add_bond_block(lines: list[str], graph: nx.Graph):
    if graph.number_of_edges() == 0:
        return

    _add_v30_line(lines, "BEGIN BOND")

    for index, edge in enumerate(graph.edges(data=True), start=1):
        node_index1, node_index2, attrs = edge
        bond_type = attrs.get("bond_type", 1)

        _add_v30_line(lines, f"{index} {bond_type} {node_index1 + 1} {node_index2 + 1}")

    _add_v30_line(lines, "END BOND")
