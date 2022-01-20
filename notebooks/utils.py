from rdkit import Chem
from rdkit.Chem import RenumberAtoms, CanonicalRankAtoms, MolFromMolFile
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from tabulate import tabulate
from collections import deque
from operator import gt, lt, eq
import random
import re
import matplotlib.pyplot as plt
import networkx as nx


def rdkit_to_nx(m):
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
    nx.set_node_attributes(nx_graph, False, "explored")
    return nx_graph

def _draw_networkx_graph(m, ax, highlight):
  nx_graph = rdkit_to_nx(m) if type(m) != nx.Graph else m
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

def write_string_representation(m):
    atom_tuples = [tuple(sorted((a.GetIdx(), n.GetIdx()))) for a in m.GetAtoms() for n in a.GetNeighbors()]
    unique_atom_tuples = sorted(set((atom_tuples)))
    return f"{CalcMolFormula(m)}:{''.join([str(t) for t in unique_atom_tuples])}"

def build_molecule_from_string_representation(string_representation):
    sum_formula, atom_tuples = string_representation.split(":")
    atom_tuples = [list(map(int, list(map(str.strip, t.split(","))))) for t in atom_tuples[1:-1].split(")(")]
    n_atoms = sum([int(s) for s in re.findall(r'-?\d+', sum_formula)])
    m = Chem.Mol()
    em = Chem.RWMol(m)
    for _ in range(n_atoms):
        em.AddAtom(Chem.Atom(1))
    for t in atom_tuples:
        em.AddBond(*t, Chem.BondType.SINGLE)
    return em.GetMol()

def permute_molecule(m, random_seed=42):
    permuted_indices = list(range(m.GetNumAtoms()))
    random.seed(random_seed)
    random.shuffle(permuted_indices)
    return RenumberAtoms(m, permuted_indices)

def print_molecule(m, caption=""):
    print(caption)
    table = []
    for atom in m.GetAtoms():
        idx = atom.GetIdx()
        num = atom.GetAtomicNum()
        partition = atom.GetIntProp("partition") if atom.HasProp("partition") else 0
        neighbors = [(n.GetIdx(),
                      n.GetAtomicNum(),
                      n.GetIntProp("partition") if n.HasProp("partition") else 0) for n in atom.GetNeighbors()]
        table.append([idx, num, partition, neighbors])
    print(tabulate(table, tablefmt="fancy_grid",
                   headers=["index", "atomic number", "partition", "neighbors (index, atomic number, partition)"]))





def _sort_molecule_by_property(m, sorting_property):
    prop_with_idcs = [(j, i) for i, j in enumerate(sorting_property)] # [(0, 0), (2, 1), (1, 2)]
    sorted_prop, idcs_sorted_by_prop = zip(*sorted(prop_with_idcs)) # (0, 1, 2), (0, 2, 1)
    return RenumberAtoms(m, idcs_sorted_by_prop)

def _atomic_number_sequence(atom):
    atomic_num = atom.GetAtomicNum()
    atomic_nums_neighbors = sorted([n.GetAtomicNum() for n in atom.GetNeighbors()], reverse=True)
    return [atomic_num] + atomic_nums_neighbors

def partition_molecule_by_atomic_numbers(m):
    '''Mutates `m`'''
    for atom in m.GetAtoms():
        atom.SetIntProp("partition", 0) # initialize partitions
    current_partition = 0
    for i in range(m.GetNumAtoms() - 1):
        j = i + 1
        atomic_numbers_i = _atomic_number_sequence(m.GetAtomWithIdx(i))
        atomic_numbers_j = _atomic_number_sequence(m.GetAtomWithIdx(j))
        if (atomic_numbers_i != atomic_numbers_j): current_partition += 1
        m.GetAtomWithIdx(j).SetIntProp("partition", current_partition)
    return m

def sort_molecule_by_canonical_indices(m):
    '''https://gist.github.com/ptosco/36574d7f025a932bc1b8db221903a8d2'''
    return _sort_molecule_by_property(m, CanonicalRankAtoms(m))

def sort_molecule_by_atomic_numbers(m):
    '''Sort atoms lexicographically by atomic number. Also, within partitions
    where all atoms have the same atomic number, sort atoms lexicographically
    by the decreasing atomic numbers of their neighbors.
    '''
    atomic_numbers = [_atomic_number_sequence(atom) for atom in m.GetAtoms()]
    return _sort_molecule_by_property(m, atomic_numbers)

def _partition_sequence(atom):
    partition = atom.GetIntProp("partition")
    partitions_neighbors = sorted([n.GetIntProp("partition") for n in atom.GetNeighbors()], reverse=True)
    return [partition] + partitions_neighbors

def _sort_molecule_by_partitions(m):
    partitions = [_partition_sequence(atom) for atom in m.GetAtoms()]
    return _sort_molecule_by_property(m, partitions)

def partition_molecule_recursively(m, show_steps=False):
    m_sorted = _sort_molecule_by_partitions(m)
    if show_steps:
        print_molecule(m_sorted, "refined partitions")
    current_partitions = [a.GetIntProp("partition") for a in m.GetAtoms()]
    updated_partitions = [0]
    for i in range(m_sorted.GetNumAtoms() - 1):
        j = i + 1
        partitions_i = _partition_sequence(m_sorted.GetAtomWithIdx(i))
        partitions_j = _partition_sequence(m_sorted.GetAtomWithIdx(j))
        current_partition = updated_partitions[-1]
        if (partitions_i != partitions_j):
            current_partition += 1
        updated_partitions.append(current_partition)
    if current_partitions == updated_partitions:
        return m_sorted
    for i, atom in enumerate(m_sorted.GetAtoms()):
        atom.SetIntProp("partition", updated_partitions[i])
    return partition_molecule_recursively(m_sorted, show_steps=show_steps)

def bfs_molecule(m, root_idx):
    """Breadth-first search over atoms.
    Note that NetworkX provides the same algorithm in `dfs_edges()`.
    This (re-)implementation allows for controlling the branching behavior
    during the molecule traversal.
    m: NetworkX graph.
    root_idx: atom at which to start traversal.
    """
    m.nodes[root_idx]["explored"] = True
    atom_queue = deque([root_idx])
    while atom_queue:
        a = atom_queue.popleft()
        for n in m.neighbors(a):
            if m.nodes[n]["explored"]:
                continue
            yield (a, n)
            m.nodes[n]["explored"] = True
            atom_queue.append(n)

def dfs_molecule(m, root_idx):
    """Depth-first search over atoms.
    Note that NetworkX provides the same algorithm in `bfs_edges()`.
    This (re-)implementation allows for controlling the branching behavior
    during the molecule traversal.
    m: NetworkX graph.
    root_idx: atom at which to start traversal.
    """
    m.nodes[root_idx]["explored"] = True
    for n_idx in m.neighbors(root_idx):
        if m.nodes[n_idx]["explored"]:
            continue
        yield (root_idx, n_idx)
        yield from dfs_molecule(m, n_idx)

def edge_dfs_molecule(m, root_idx):
    """Depth-first search over edges.
    Note that NetworkX provides the same algorithm in `edge_dfs ()`.
    This (re-)implementation allows for controlling the branching behavior
    during the molecule traversal.
    m: NetworkX graph.
    root_idx: atom at which to start traversal.
    """
    visited_edges = set()
    visited_nodes = set()
    edges = {}

    nodes = list(m.nbunch_iter(root_idx))
    for start_node in nodes:
        stack = [start_node]
        while stack:
            current_node = stack[-1]
            if current_node not in visited_nodes:
                edges[current_node] = iter(m.edges(current_node))
                visited_nodes.add(current_node)

            try:
                edge = next(edges[current_node])
            except StopIteration:
                # No more edges from the current node.
                stack.pop()
            else:
                edgeid = (frozenset(edge[:2]),) + edge[2:]
                if edgeid not in visited_edges:
                    visited_edges.add(edgeid)
                    # Mark the traversed "to" node as to-be-explored.
                    stack.append(edge[1])
                    yield edge


def test_invariance(m, canonicalization_steps=[]):
    """`canonicalization_steps`: list of functions that transform initial
    molecule into canonicalized molecule. Functions are applied in the order in
    which they appear in the list."""
    m_permuted = permute_molecule(m)

    for func in canonicalization_steps:
      m = func(m)
      m_permuted = func(m_permuted)

    string_original = write_string_representation(m)
    string_permuted = write_string_representation(m_permuted)

    try:
        assert string_original == string_permuted, f"{string_original}\ndoesn't match\n{string_permuted}."
        print(f"{string_original}\nmatches\n{string_permuted}.")
        return True
    except AssertionError as e:
        print(e)
        return False

def load_molfile(path):
    m = MolFromMolFile(str(path), removeHs=False)
    try:
        # Unfortunately, RDKit exceptions cannot be caught. However, MolFromMolFile returns None if validation fails.
        n_atoms = m.GetNumAtoms() # Raises AttributeError in case m_original is None.
        if n_atoms < 1:
            raise Exception(f"Molecule has < 1 ({n_atoms}) atoms")
        return m
    except Exception as e:
#         print(f"Failed loading {str(molfile)}: {e}.")
        return None

def partition_molecule(m):
    m_sorted_by_atomic_numbers = sort_molecule_by_atomic_numbers(m)
    m_partitioned_by_atomic_numbers = partition_molecule_by_atomic_numbers(m_sorted_by_atomic_numbers)
    return partition_molecule_recursively(m_partitioned_by_atomic_numbers, show_steps=False)

def canonicalize_molecule(m, root_idx=0):
    m_partitioned = partition_molecule(m)
    m_partitioned = rdkit_to_nx(m_partitioned) if type(m) != nx.Graph else m
    canonical_idcs = traverse_molecule(m_partitioned, root_idx)
    return nx.relabel_nodes(m_partitioned, canonical_idcs, copy=True)

def create_partition_lut(m):
    """Look-up-table of atom indices in partitions."""
    partitions = set(sorted([v for k, v in m.nodes.data("partition")]))
    partition_lut = {p:set() for p in partitions}
    for a in m:
        partition_lut[m.nodes[a]["partition"]].add(a)
    partition_lut.update((k, sorted(list(v), reverse=True)) for k, v in partition_lut.items())
    return partition_lut

def traverse_molecule(m, root_idx, traversal_priorities=[lt, gt, eq], show_traversal_order=False):
    partitions = m.nodes.data("partition")
    lut = create_partition_lut(m)
    atom_stack = [root_idx]
    canonical_idcs = {}

    while atom_stack:
        a = atom_stack.pop()
        if m.nodes[a]["explored"]:
            continue
        a_canon = lut[partitions[a]].pop()
        canonical_idcs[a] = a_canon
        if show_traversal_order:
            print(f"Current atom index: {a}.\tRe-labeling to {a_canon}.")
        neighbors = list(m.neighbors(a))
        neighbor_traversal_order = []
        for priority in traversal_priorities:
            neighbors_this_priority = [n for n in neighbors
                                       if priority(partitions[a], partitions[n])]
            neighbor_traversal_order.extend(sorted(neighbors_this_priority))

        m.nodes[a]["explored"] = True
        for n in neighbor_traversal_order:
            atom_stack.insert(0, n)

    return canonical_idcs
