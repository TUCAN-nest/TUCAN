# (c) CC BY-SA | Ulrich Schatzschneider | Universität Würzburg | NFDI4Chem | v1.4 | 06.06.2021

def read_molfile(filename)
  raise "#{filename} doesn't exist." unless File.exist?(filename)

  molfile = File.read(filename) # reads entire file and closes it
  molfile_lines = molfile.split("\n")

  puts "\nPrinting molfile: #{filename}. First 4 lines contain header."
  puts molfile

  molfile_lines
end

def create_molecule_array(molfile_lines, periodic_table_elements)
  # number of atoms in file is 1st number of 4th line (with index "3" as counting starts at zero)
  # number of bonds is 2nd number of 4th line (NOT: length of file minus header length (4 lines) minus number of atom definition lines { minus one (as final line is "M END") } )
  # (This is not true, since there can be additional lines with definitions starting with letter "M", thus rather use second number which is the number of bonds)
  atom_count, bond_count = molfile_lines[3].scan(/\d+/).map { |n| n.to_i }

  atomic_masses = []
  (4..atom_count + 3).each do |atom_index|
    atom = molfile_lines[atom_index].split(' ')
    atomic_masses.push(periodic_table_elements.index(atom[3]))
  end

  molecule_graph = []
  (0..atom_count - 1).each do |atom_index|
    molecule_graph.push([atom_index, []])
  end

  # now read the remaining lines containing the bond definitions in the sequence atom1 atom2 bond_order ... unknown/unused ... (can be ignored)
  (0..bond_count - 1).each do |bond_index|
    vertex1, vertex2, * = molfile_lines[bond_index + 4 + atom_count].split(' ').map { |i| i.to_i - 1 }
    vertex1, vertex2 = vertex2, vertex1 if vertex1 > vertex2 # make sure first atom always has lower (not: higher?) index
    molecule_graph[vertex1][1].push(vertex2)    # need to push twice, to the first atom of a bond
    molecule_graph[vertex2][1].push(vertex1)    # and then to the second atom of the bond
  end

  molecule_graph.map! do |atom|
    atom, edges = atom
    atom = [atom, atomic_masses[atom]]
    edges.map! { |edge| [edge, atomic_masses[edge]]}
    [atom, edges]
  end

  molecule_graph
end

def canonicalize_molecule(molecule)
  puts "\nInitial molecule:"
  molecule.each { |atom| puts atom.inspect }

  sorted_molecule = sort_elements_by_atomic_mass(molecule)
  puts "\nMolecule with elements sorted by atomic mass (increasing):"
  sorted_molecule.each { |atom| puts atom.inspect }

  sorted_molecule = update_molecule_indices(sorted_molecule, update_edges=true)
  puts "\nMolecule with updated indices after sorting:"
  sorted_molecule.each { |atom| puts atom.inspect }

  sorted_molecule = sort_edges_by_atomic_mass(sorted_molecule)
  puts "\nMolecule with edges sorted by atomic mass (decreasing) per element:"
  sorted_molecule.each { |atom| puts atom.inspect }

  sort_elements_by_number_of_edges(sorted_molecule)
  puts "\nMolecule with elements of same kind sorted by number of edges (increasing):"
  sorted_molecule.each { |atom| puts atom.inspect }

  sorted_molecule = update_molecule_indices(sorted_molecule, update_edges=true)
  puts "\nMolecule with updated indices after sorting:"
  sorted_molecule.each { |atom| puts atom.inspect }

  sort_elements_by_index_of_edges(sorted_molecule)
  puts "\nMolecule with elements of same kind and same number of edges sorted by atomic mass of edges (increasing):"
  sorted_molecule.each { |atom| puts atom.inspect }

  sorted_molecule = update_molecule_indices(sorted_molecule, update_edges=true)
  puts "\nMolecule with updated indices after sorting:"
  sorted_molecule.each { |atom| puts atom.inspect }

  sorted_molecule
end

def write_ninchi_string(molecule, periodic_table_elements)
  sum_formula = write_sum_formula_string(molecule, periodic_table_elements)
  serialized_molecule = serialize_molecule(molecule)
  "nInChI=1S/#{sum_formula}/c#{serialized_molecule}"
end

def write_dot_file(molecule, periodic_table_elements, periodic_table_colors)
  dotfile = "graph test\n{\n  bgcolor=grey\n"
  molecule.each_with_index do |atom, i|
    symbol = periodic_table_elements[atom[0][1]]
    color = periodic_table_colors.fetch(symbol, 'lightgrey')
    dotfile += "  #{i} [label=\"#{symbol} #{i}\" color=#{color},style=filled,shape=circle,fontname=Calibri];\n"
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
    element, edges = atom
    edges.each do |edge|
      graph.push([element[0], edge[0]].sort)
    end
  end
  graph.uniq.sort
end

def write_sum_formula_string(molecule, periodic_table_elements)
  # Write sum formula in the order C > H > all other elements in alphabetic order.
  element_counts = compute_element_counts(molecule, periodic_table_elements)
  element_counts.transform_values! { |v| v > 1 ? v : '' } # remove 1s since counts of 1 are implicit in sum formula
  sum_formula_string = ''
  sum_formula_string += "C#{element_counts['C']}" if element_counts.key?('C')
  sum_formula_string += "H#{element_counts['H']}" if element_counts.key?('H')
  element_counts.sort.to_h.each do |element, count|
    sum_formula_string += "#{element}#{count}" unless %w[C H].include?(element)
  end
  sum_formula_string
end

def compute_element_counts(molecule, periodic_table_elements)
  # Compute hash table mapping element symbols to stoichiometric counts.
  unique_elements = molecule.map { |atom| atom[0][1] }.uniq
  initial_counts = Array.new(unique_elements.length, 0)
  element_counts = unique_elements.zip(initial_counts).to_h
  molecule.each { |atom| element_counts[atom[0][1]] += 1 }
  element_counts.transform_keys! { |k| periodic_table_elements[k] } # change atomic mass to element symbol
end

def sort_edges_by_atomic_mass(molecule)
  # Sort edge numbers of each atom left to right from large to small.
  sorted_molecule = []
  molecule.each do |atom|
    element, edges = atom
    sorted_edges = edges.sort { |a, b| b[1] <=> a[1] }
    sorted_molecule.push([element, sorted_edges])
  end
  sorted_molecule
end

def sort_elements_by_index_of_edges(molecule)
  # Mutates `molecule`.
  molecule.sort! do |atom_a, atom_b|
    mass_a, mass_b = atom_a[0][1], atom_b[0][1]
    mass_edges_a = atom_a[1].map { |atom| atom[1] }
    mass_edges_b = atom_b[1].map { |atom| atom[1] }
    case
    when (mass_a == mass_b) && # A and B are the same element ...
         (mass_edges_a.length == mass_edges_b.length) && # with the same number of edges ...
         (mass_edges_a <=> mass_edges_b) == 1 # and the edge indices of A are larger than the ones of B.
      1
    when (mass_a == mass_b) && # A and B are the same element ...
         (mass_edges_a.length == mass_edges_b.length) && # with the same number of edges ...
         (mass_edges_a <=> mass_edges_b) == -1 # and the edge indices of A are smaller than the ones of B.
      -1
    else
      0
    end
  end
end

def sort_elements_by_number_of_edges(molecule)
  # Mutates `molecule`.
  molecule.sort! do |atom_a, atom_b|
    case
    when (atom_a[0][1] == atom_b[0][1]) && # A and B are the same element ...
         (atom_a[1].length > atom_b[1].length) # and A has more edges than B.
      1
    when (atom_a[0][1] == atom_b[0][1]) && # A and B are the same element ...
         (atom_a[1].length < atom_b[1].length) # and A has less edges than B.
      -1
    else
      0
    end
  end
end

def sort_elements_by_atomic_mass(molecule)
  molecule.sort { |a, b| a[0][1] <=> b[0][1] }
end

def update_molecule_indices(molecule, update_edges=false)
  index_updates = compute_index_updates(molecule)
  updated_molecule = []
  molecule.each do |atom|
    element, edges = atom
    updated_element = update_element_index(element, index_updates)
    updated_edges = update_edges ? update_edge_indices(edges, index_updates) : edges
    updated_molecule.push([updated_element, updated_edges])
  end
  updated_molecule
end

def update_element_index(element, index_updates)
  element[0] = index_updates[element[0]] if index_updates.key?(element[0])
  element
end

def update_edge_indices(edges, index_updates)
  edges.map do |edge|
    edge[0] = index_updates[edge[0]] if index_updates.key?(edge[0])
    edge
  end
end

def compute_index_updates(molecule)
  index_updates = {}
  molecule.each_with_index do |atom, i|
    element, * = atom
    index_updates[element[0]] = i if element[0] != i
  end
  index_updates
end
