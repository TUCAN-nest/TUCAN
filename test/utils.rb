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

  def ninchi_string(permute_molfile: false)
    atom_block, edge_block = molfile_data(permute_molfile)
    adjacency_matrix, node_features_matrix = create_adjacency_matrix(atom_block, edge_block, PeriodicTable::ELEMENTS)
    adjacency_matrix, node_features_matrix = sort_adjacency_matrix(adjacency_matrix, node_features_matrix)
    write_ninchi_string(adjacency_matrix, node_features_matrix, PeriodicTable::ELEMENTS)
  end

  def molfile_data(permute)
    atom_block, edge_block = read_molfile(@path)
    atom_block, edge_block = permute_molfile_data(atom_block, edge_block) if permute
    [atom_block, edge_block]
  end

  def permute_molfile_data(atom_block, edge_block)
    puts 'permuting'
    [atom_block, edge_block]
  end
end
