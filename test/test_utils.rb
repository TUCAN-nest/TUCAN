require 'minitest/autorun'
require './test/utils'

class UtilsTests < Minitest::Test
  # Test if permutation works correctly.
  # Permuting ammonia with random seed 181 results in the following changes.
  # Atom block:
  # original: ["M V30 1 H 0.0 0.0 0.0000 0", "M V30 2 H 0.0 0.0 0.0000 0", "M V30 3 H 0.0 0.0 0.0000 0", "M V30 4 N 0.0 0.0 0.0000 0"]
  # permuted: ["M V30 1 H 0.0 0.0 0.0000 0", "M V30 2 H 0.0 0.0 0.0000 0", "M V30 3 N 0.0 0.0 0.0000 0", "M V30 4 H 0.0 0.0 0.0000 0"]
  # Edge block
  # original: ["M  V30 1 1 4 1", "M  V30 2 1 2 4", "M  V30 3 1 4 3"]
  # permuted: ["M V30 1 1 3 1", "M V30 2 1 4 3", "M V30 3 1 3 2"]

  def setup
    @ammonia = Molecule.new('test/testfiles/ammonia/ammonia.mol')
  end

  def test_permutation_shuffles_atom_order
    # Assert that the atom block of ammonia is permuted.
    atom_block, = @ammonia.molfile_data(false, 181)
    permuted_atom_block, = @ammonia.molfile_data(true, 181)
    atoms = atom_block.map { |atom| atom.split(' ')[3] }
    permuted_atoms = permuted_atom_block.map { |atom| atom.split(' ')[3] }
    assert_equal(%w[H H H N], atoms)
    assert_equal(%w[H H N H], permuted_atoms)
  end

  def test_permutation_preserves_atom_indices
    # Assert that the indices of the atom block of ammonia are preserved after permutation.
    # That is, only the order of the atom symbols should be permuted, not the order of the indices.
    atom_block, = @ammonia.molfile_data(false, 181)
    permuted_atom_block, = @ammonia.molfile_data(true, 181)
    atom_indices = atom_block.map { |atom| atom.split(' ')[2] }
    permuted_atom_indices = permuted_atom_block.map { |atom| atom.split(' ')[2] }
    assert_equal(atom_indices, permuted_atom_indices)
  end

  def test_permutation_updates_bonds
    # Assert that the bond vertices in the bond block of ammonia are updated after permutation.
    _, edge_block = @ammonia.molfile_data(false, 181)
    _, permuted_edge_block = @ammonia.molfile_data(true, 181)
    edges = edge_block.map { |edge| edge.split(' ')[4..5] }
    permuted_edges = permuted_edge_block.map { |edge| edge.split(' ')[4..5] }
    assert_equal([['4', '1'], ['2', '4'], ['4', '3']], edges) # all H's connected to N
    assert_equal([['3', '1'], ['4', '3'], ['3', '2']], permuted_edges) # unchanged, all H's connected to N
  end

  def test_permutation_preserves_bond_indices
    # Assert that the indices of the bond block of ammonia are preserved after permutation.
    # That is, only the bond vertices should be updated according to the atom permutations,
    # while leaving the rest of the bond block unaffected.
    _, edge_block = @ammonia.molfile_data(false, 181)
    _, permuted_edge_block = @ammonia.molfile_data(true, 181)
    edge_indices = edge_block.map { |edge| edge.split(' ')[2] }
    permuted_edge_indices = permuted_edge_block.map { |edge| edge.split(' ')[2] }
    assert_equal(edge_indices, permuted_edge_indices)
  end
end
