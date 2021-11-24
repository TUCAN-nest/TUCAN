require './inchi'
require './periodic_table'

class Molecule
  attr_reader :reference_ninchi_string, :name

  def initialize(path, ninchi_string = nil)
    @path = path
    @name = File.basename(path, '.mol')
    @reference_ninchi_string = ninchi_string
  end

  def ninchi_string(permute_molfile: false, random_seed: Random.new_seed)
    molfile_lines = molfile_lines(permute_molfile, random_seed)
    adjacency_matrix, node_features_matrix, distance_matrix, molfile_header, molfile_footer = initialize_matrix(molfile_lines, PeriodicTable::ELEMENTS)
    adjacency_matrix, node_features_matrix, distance_matrix = sort_adjacency_matrix(adjacency_matrix, node_features_matrix, distance_matrix)
    write_ninchi_string(adjacency_matrix, node_features_matrix, PeriodicTable::ELEMENTS)
  end

  def molfile_lines(permute, random_seed)
    molfile_lines = read_molfile(@path)
    molfile_lines = permute_molfile(molfile_lines, random_seed) if permute
    molfile_lines
  end

  def permute_molfile(molfile_lines, random_seed)
    atom_count = molfile_lines[5].split(' ')[3].to_i
    edge_count = molfile_lines[5].split(' ')[4].to_i
    molfile_lines = permute_atoms(molfile_lines, atom_count, random_seed)
    permuted_indices = (1..atom_count).to_a.shuffle(random: Random.new(random_seed))
    permute_edges(molfile_lines, atom_count, edge_count, permuted_indices)
  end

  def permute_atoms(molfile_lines, atom_count, random_seed)
    start_atom_block = 7
    end_atom_block = start_atom_block + atom_count
    atom_block = (start_atom_block...end_atom_block).map { |i| molfile_lines[i] }
    atom_block.shuffle!(random: Random.new(random_seed))
    atom_block = update_atom_ids(atom_block)
    molfile_lines[start_atom_block...end_atom_block] = atom_block
    molfile_lines
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

  def permute_edges(molfile_lines, atom_count, edge_count, permuted_indices)
    start_edge_block = 9 + atom_count
    end_edge_block = start_edge_block + edge_count
    edge_block = (start_edge_block...end_edge_block).map { |i| molfile_lines[i] }
    edge_block = update_edge_ids(edge_block, permuted_indices)
    molfile_lines[start_edge_block...end_edge_block] = edge_block
    molfile_lines
  end

  def update_edge_ids(edge_block, permuted_indices)
    updated_edge_block = []
    edge_block.each do |edge|
      updated_edge = update_edge(edge, permuted_indices)
      updated_edge_block.push(updated_edge)
    end
    updated_edge_block
  end

  def update_edge(edge, permuted_indices)
    updated_edge = edge.split(' ')
    updated_edge[4] = permuted_indices.index(updated_edge[4].to_i) + 1
    updated_edge[5] = permuted_indices.index(updated_edge[5].to_i) + 1
    updated_edge.join(' ')
  end
end
