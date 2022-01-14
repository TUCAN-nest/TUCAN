from rdkit import Chem
from rdkit.Chem import RenumberAtoms, CanonicalRankAtoms, MolFromMolFile
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from tabulate import tabulate
from collections import deque
import random
import re
import matplotlib.pyplot as plt
import networkx as nx


def _draw_networkx_graph(m, ax, highlight):
  edges = [(b.GetBeginAtom().GetIdx(), b.GetEndAtom().GetIdx()) for b in m.GetBonds()]
  nx_graph = nx.from_edgelist(edges)
  attribute = [a.GetAtomicNum() for a in m.GetAtoms()]
  if highlight == "partition":
    attribute = [a.GetIntProp("partition") if a.HasProp("partition") else 0 for a in m.GetAtoms()]
  nx.set_node_attributes(nx_graph, {i:p for i, p in enumerate(attribute)}, highlight)
  highlight_colors = list(nx.get_node_attributes(nx_graph, highlight).values())
  node_size = 1 / nx_graph.order() * 10000
  nx.draw_kamada_kawai(nx_graph, node_color=highlight_colors, node_size=node_size,
                       cmap="rainbow", alpha=.5, with_labels=True, font_weight="heavy", ax=ax)

def draw_molecules(m_list, caption_list, highlight="atomic_number"):
  """`highlight`: color atoms by "atomic_number" (default) or "partition"."""
  if highlight not in ["atomic_number", "partition"]:
    print("Please select one of {'partition', 'atomic_number'} for `highlight`.")
    return
  n_molecules = len(m_list)
  fig = plt.figure(figsize=(n_molecules * 6, 6))
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

def permute_molecule(m):
    permuted_indices = list(range(m.GetNumAtoms()))
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

def bfs_molecule(m):
    for a in m.GetAtoms():
        a.SetBoolProp("explored", False)
    root_atom = m.GetAtomWithIdx(0)
    root_atom.SetBoolProp("explored", True)
    atom_queue = deque([root_atom])
    while atom_queue:
        atom = atom_queue.popleft()
        for n in atom.GetNeighbors():
            if n.GetBoolProp("explored"):
                continue
            n_neighbors = n.GetNeighbors()
            if len([nn.GetIntProp("partition") for nn in n_neighbors]) != len(n_neighbors):
                continue
            print(atom.GetIdx(), n.GetIdx())
            n.SetBoolProp("explored", True)
            atom_queue.append(n)

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


