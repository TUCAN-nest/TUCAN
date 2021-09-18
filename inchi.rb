# (c) CC BY-SA | Ulrich Schatzschneider | Universität Würzburg | NFDI4Chem | v1.4 | 06.06.2021

require './periodic_table'

class Atom
  include PeriodicTable
  attr_reader :mass, :symbol, :color
  attr_accessor :id
  attr_writer :edges

  def initialize(id, symbol)
    @id = id
    @symbol = symbol
    @mass = PeriodicTable::ELEMENTS.index(@symbol)
    @color = PeriodicTable::ELEMENT_COLORS[@symbol]
    @edges = []
  end

  def edges
    @edges.sort!.reverse! # return in descending order
  end
end

module Inchi

  def read_molfile(filename)
    if(filename.nil?)
      print "\nPlease provide a filename.\n"
      exit(false)
    end
    unless File.exist?(filename)
      print "\n#{filename} doesn't exist.\n"
      exit(false)
    end
    
    molfile = File.read(filename) # reads entire file and closes it
    molfile.split("\n")
  end

  def create_molecule_array(molfile_lines)
    molecule = []
    atom_count, edge_count = molfile_lines[3].scan(/\d+/).map(&:to_i) # on 4th  line, 1st number is number of atoms, 2nd number is number of bonds.
    (4..atom_count + 3).each_with_index do |atom_index, i|
      molecule.push(Atom.new(i, molfile_lines[atom_index].split(' ')[3]))
    end
    (0..edge_count - 1).each do |edge_index|
      vertex1, vertex2 = parse_edge(molfile_lines[edge_index + 4 + atom_count])
      molecule[vertex1].edges.push(vertex2)    # add to the first atom of a bond
      molecule[vertex2].edges.push(vertex1)    # and to the second atom of the bond
    end
    molecule
  end

  def canonicalize_molecule(molecule, filename)
    filename = File.basename(filename, '.mol')

    print_molecule(molecule,
      "\nInitial data structure of #{filename}:")

    sorted_molecule = sort_elements_by_atomic_mass(molecule)
    print_molecule(sorted_molecule,
      "\n#{filename} with atoms sorted by atomic mass (increasing):")

    sorted_molecule = sort_elements_by_number_of_edges(sorted_molecule)
    print_molecule(sorted_molecule,
      "\n#{filename} with atoms of same kind sorted by number of edges (increasing):")

    sorted_molecule = update_molecule_indices(sorted_molecule)
    print_molecule(sorted_molecule,
      "\n#{filename} with updated indices after sorting:")

    *previous_molecule_states, sorted_molecule = sort_elements_by_index_of_edges(sorted_molecule)
    print_molecule(sorted_molecule,
      "\n#{filename} with atoms of same kind and same number of edges sorted by indices of edges (increasing):")

    inspect_molecule_states(previous_molecule_states, sorted_molecule, filename)

    sorted_molecule
  end

  def write_ninchi_string(molecule)
    sum_formula = write_sum_formula_string(molecule)
    serialized_molecule = serialize_molecule(molecule)
    "nInChI=1S/#{sum_formula}/c#{serialized_molecule}"
  end

  def write_dot_file(molecule, filename)
    filename = File.basename(filename, '.mol')
    dotfile = "graph #{filename}\n{\n  bgcolor=grey\n"
    molecule.each_with_index do |atom, i|
      dotfile += "  #{i} [label=\"#{atom.symbol} #{i}\" color=#{atom.color},style=filled,shape=circle,fontname=Calibri];\n"
    end
    graph = compute_graph(molecule)
    graph.each { |line| dotfile += "  #{line[0]} -- #{line[1]} [color=black,style=bold];\n" if line[0] != line[1] }
    dotfile += "}\n"
  end

  # Helper methods ###################################################################################
  private

  def serialize_molecule(molecule)
    graph = compute_graph(molecule)
    inchi_string = ''
    graph.each do |line|
      inchi_string += "(#{line[0]}-#{line[1]})" if line[0] != line[1]
    end
    inchi_string
  end

  def compute_graph(molecule)
    graph = []
    molecule.each do |atom|
      atom.edges.each do |edge|
        graph.push([atom.id, edge].sort)
      end
    end
    graph.uniq.sort
  end

  def write_sum_formula_string(molecule)
    # Write sum formula in the order C > H > all other elements in alphabetic order.
    element_counts = compute_element_counts(molecule)
    element_counts.transform_values! { |v| v > 1 ? v : '' } # remove 1s since counts of 1 are implicit in sum formula
    sum_formula_string = ''
    sum_formula_string += "C#{element_counts['C']}" if element_counts.key?('C')
    sum_formula_string += "H#{element_counts['H']}" if element_counts.key?('H')
    element_counts.sort.to_h.each do |element, count|
      sum_formula_string += "#{element}#{count}" unless %w[C H].include?(element)
    end
    sum_formula_string
  end

  def compute_element_counts(molecule)
    # Compute hash table mapping element symbols to stoichiometric counts.
    unique_elements = molecule.map { |atom| atom.symbol }.uniq
    initial_counts = Array.new(unique_elements.length, 0)
    element_counts = unique_elements.zip(initial_counts).to_h
    molecule.each { |atom| element_counts[atom.symbol] += 1 }
    element_counts
  end

  def sort_elements_by_index_of_edges(molecule)
    # Cannot use built-in sort since indices have to be updated after every swap,
    # rather than once after sorting is done (like with the other sorting steps).
    # This is because we sort by indices.
    n_iterations = molecule.size - 2
    previous_molecule_states = [Marshal.load(Marshal.dump(molecule))]
    sorted = false
    until sorted

      for i in 0..n_iterations
        atom_a = molecule[i]
        atom_b = molecule[i + 1]

        # Swap A and B (i.e., bubble up A) if ...
        if (atom_a.mass == atom_b.mass) && # A and B are the same element ...
          (atom_a.edges.length == atom_b.edges.length) && # with the same number of edges ...
          (atom_a.edges <=> atom_b.edges) == 1 # and A is connected to larger indices than B.
          # Spaceship operator (<=>) compares arrays pairwise element-by-element.
          # I.e., first compare the two elements at index 0, etc.. Result is determined by first unequal element pair.
          # The operator returns 1 if A > B, -1 if A < B, and 0 if A == B.
          molecule[i], molecule[i + 1] = molecule[i + 1], molecule[i]
          molecule = update_molecule_indices(molecule)
        end
      end
      sorted = previous_molecule_states.map { |state| state.map(&:edges) }.include?(molecule.map(&:edges))

      previous_molecule_states.push(Marshal.load(Marshal.dump(molecule)))
    end
    previous_molecule_states
  end

  def sort_elements_by_number_of_edges(molecule)
    # Note that `each_with_index` returns position of atom in molecule array,
    # not atom ID.
    molecule.each_with_index.sort { |(atom_a, idx_a), (atom_b, idx_b)|
      order = 0 # assume unequal masses and hence no swap; note that sort is not stable if order = 0
      order = atom_a.edges.length <=> atom_b.edges.length if atom_a.mass == atom_b.mass # sort by vertex degree in case of equal masses
      order = idx_a <=> idx_b if order.zero? # sort by index in case of a) unequal masses or b) equal masses and equal vertex degree; this ensures stable sort
      order
    }.map(&:first) # discard index (second element)
  end

  def sort_elements_by_atomic_mass(molecule)
    # Note that `each_with_index` returns position of atom in molecule array,
    # not atom ID.
    molecule.each_with_index.sort { |(atom_a, idx_a), (atom_b, idx_b)|
      order = atom_a.mass <=> atom_b.mass
      order = idx_a <=> idx_b if order.zero? # sort by index in case of equal masses (i.e., order = 0); this ensures stable sort
      order
    }.map(&:first) # discard index (second element)
  end

  def update_molecule_indices(molecule, random_indices=false)
    index_updates = compute_index_updates(molecule, random_indices)
    updated_molecule = Marshal.load(Marshal.dump(molecule))
    updated_molecule.each do |atom|
      atom.id = index_updates[atom.id]
      atom.edges.map! { |edge| index_updates[edge] }
    end
    updated_molecule
  end

  def compute_index_updates(molecule, random_indices)
    current_indices = molecule.map(&:id)
    updated_indices = (0..molecule.length - 1).to_a
    updated_indices.shuffle! if random_indices
    current_indices.zip(updated_indices).to_h
  end

  def parse_edge(molfile_line)
    vertex1, vertex2, * = molfile_line.split(' ').map { |i| i.to_i - 1 }
    vertex1, vertex2 = vertex2, vertex1 if vertex1 > vertex2    # make sure first atom always has lower (not: higher?) index
    [vertex1, vertex2]
  end

  def print_molecule(molecule, caption)
    puts caption
    puts "\nindex\tmass\tindices of connected atoms"
    puts "-----\t----\t--------------------------"
    molecule.each { |atom| puts "#{atom.id}\t#{atom.mass + 1}\t#{atom.edges}" }
  end

  def inspect_molecule_states(previous_states, final_state, filename)
    puts "\nPrinting all molecule states of #{filename} that occured during sorting by indices of edges..."
    previous_states.each_with_index do |state, i|
      print_molecule(state,
        "\nIteration #{i} yielded the following state:")
    end
    r = previous_states.map { |state| state.map(&:edges) }.index(final_state.map(&:edges))
    print_molecule(final_state,
      "\nSorting converged in iteration #{previous_states.size} with re-occurence of state at iteration #{r}:")
  end

end
