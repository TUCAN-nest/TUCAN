require './inchi'
require './periodic_table'

class Molecule
  include Inchi
  attr_reader :reference_ninchi_string, :name

  def initialize(path, ninchi_string = nil)
    @path = path
    @name = File.basename(path, '.mol')
    @reference_ninchi_string = ninchi_string
  end

  def ninchi_string(permute_molfile: false, random_seed: Random.new_seed)
    atom_block, edge_block = molfile_data(permute_molfile, random_seed)
    adjacency_matrix, node_features_matrix = create_adjacency_matrix(atom_block, edge_block, PeriodicTable::ELEMENTS)
    adjacency_matrix, node_features_matrix = sort_adjacency_matrix(adjacency_matrix, node_features_matrix)
    write_ninchi_string(adjacency_matrix, node_features_matrix, PeriodicTable::ELEMENTS)
  end

  def molfile_data(permute, random_seed)
    atom_block, edge_block = read_molfile(@path)
    atom_block, edge_block = permute_molfile_data(atom_block, edge_block, random_seed) if permute
    [atom_block, edge_block]
  end

  def permute_molfile_data(atom_block, edge_block, random_seed)
    atom_block = permute_atoms(atom_block, random_seed)
    permuted_indices = (1..atom_block.size).to_a.shuffle(random: Random.new(random_seed))
    edge_block = permute_edges(edge_block, permuted_indices)
    [atom_block, edge_block]
  end

  def permute_atoms(atom_block, random_seed)
    atom_block.shuffle!(random: Random.new(random_seed))
    update_atom_ids(atom_block)
  end

  def update_atom_ids(atom_block)
    updated_atom_block = []
    atom_block.each.with_index do |atom, i|
      updated_atom = atom.split(' ')
      updated_atom[2] = i + 1
      updated_atom_block.push(updated_atom.join(' '))
    end
    updated_atom_block
  end

  def permute_edges(edge_block, permuted_indices)
    permuted_edge_block = []
    edge_block.each do |edge|
      updated_edge = update_edge_ids(edge, permuted_indices)
      permuted_edge_block.push(updated_edge)
    end
    permuted_edge_block
  end

  def update_edge_ids(edge, permuted_indices)
    updated_edge = edge.split(' ')
    updated_edge[4] = permuted_indices.index(updated_edge[4].to_i) + 1
    updated_edge[5] = permuted_indices.index(updated_edge[5].to_i) + 1
    updated_edge.join(' ')
  end
end
