from rdkit.Chem import RenumberAtoms, CanonicalRankAtoms, MolFromMolFile
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import random


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

def _sort_molecule_by_property(m, sorting_property):
    prop_with_idcs = [(j, i) for i, j in enumerate(sorting_property)] # [(0, 0), (2, 1), (1, 2)]
    sorted_prop, idcs_sorted_by_prop = zip(*sorted(prop_with_idcs)) # (0, 1, 2), (0, 2, 1)
    return RenumberAtoms(m, idcs_sorted_by_prop)

def _atomic_number_sequence(atom):
    atomic_num = atom.GetAtomicNum()
    atomic_nums_neighbors = sorted([n.GetAtomicNum() for n in atom.GetNeighbors()], reverse=True)
    return [atomic_num] + atomic_nums_neighbors

def write_string_representation(m):
    atom_tuples = [tuple(sorted((a.GetIdx(), n.GetIdx()))) for a in m.GetAtoms() for n in a.GetNeighbors()]
    unique_atom_tuples = sorted(set((atom_tuples)))
    return f"{CalcMolFormula(m)}:{''.join([str(t) for t in unique_atom_tuples])}"

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

def permute_molecule(m, random_seed=42):
    permuted_indices = list(range(m.GetNumAtoms()))
    random.seed(random_seed)
    random.shuffle(permuted_indices)
    return RenumberAtoms(m, permuted_indices)
