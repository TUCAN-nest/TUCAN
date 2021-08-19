# (c) CC BY-SA | Ulrich Schatzschneider | Universität Würzburg | NFDI4Chem | v1.4 | 06.06.2021

def read_molfile(filename)
  raise "#{filename} doesn't exist." unless File.exist?(filename)

  molfile = File.read(filename) # reads entire file and closes it
  molfile_lines = molfile.split("\n")

  puts "Start of file: #{filename}\n"
  puts "-----------------------------------"
  puts molfile
  puts "-----------------------------------"

  puts "Printing molfile header"
  puts "-----------------------------------"
  (0..3).each do |i| # header of a molfile is 4 lines long
    puts "#{i}: #{molfile_lines[i]}\n"
  end
  puts "-----------------------------------"

  molfile_lines
end

def create_molecule_array(molfile_lines, periodic_table_elements)

  # number of atoms in file is 1st number of 4th line (with index "3" as counting starts at zero)
  # number of bonds is 2nd number of 4th line (NOT: length of file minus header length (4 lines) minus number of atom definition lines { minus one (as final line is "M END") } )
  # (This is not true, since there can be additional lines with definitions starting with letter "M", thus rather use second number which is the number of bonds)
  atom_count, bond_count = molfile_lines[3].scan(/\d+/).map { |n| n.to_i }
  puts "Molecule has #{atom_count} atoms.\n"
  puts "Molecule has #{bond_count} bonds.\n"

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
    puts "Bond #{bond_index + 1}: #{vertex1}-#{vertex2}"
    molecule_graph[vertex1][1].push(vertex2)    # need to push twice, to the first atom of a bond
    molecule_graph[vertex2][1].push(vertex1)    # and then to the second atom of the bond
  end

  molecule_graph.map! do |atom|
    atom, connections = atom
    atom = [atom, atomic_masses[atom]]
    connections.map! { |connection| [connection, atomic_masses[connection]]}
    [atom, connections]
  end

  molecule_graph
end

# def calculate_sum_formula(molecule, periodic_table_elements)
#   # Create sum formula array in the order C > H > all other elements in alphabetic order.
#   # "atom count" is set to zero.
#   sum_formula = []
#   sum_formula.push(['C', 0], ['H', 0])
#   periodic_table_elements = periodic_table_elements - ['C'] - ['H']
#   periodic_table_elements.sort!
#   periodic_table_elements.each do |element|
#     sum_formula.push([element, 0])
#   end
#   sum_formula.each do |element|
#     molecule.each do |atom|
#       element[1] += 1 if element[0] == atom[0] # sum up
#     end
#   end

#   # generate output string
#   sum_formula_string = ''
#   sum_formula.each do |element|
#     if element[1] > 1
#       sum_formula_string = sum_formula_string + element[0].to_s + String(element[1])
#     elsif element[1] > 0
#       sum_formula_string += element[0].to_s # if there is only one atom of a particular element in the molecule, only add the element symbol without stoichiometric count
#     end
#   end
#   raise 'ERROR - This *.mol file does not contain any molecular structure' if sum_formula_string == ''

#   puts "The sum formula of this molecule is: #{sum_formula_string}\n"
#   print sum_formula, "\n\n"
#   sum_formula_string
# end
def compute_element_counts(molecule, periodic_table_elements)
  unique_elements = molecule.map { |atom| atom[0][1] }.uniq
  initial_counts = Array.new(unique_elements.length, 0)
  element_counts = unique_elements.zip(initial_counts).to_h
  molecule.each { |atom| element_counts[atom[0][1]] += 1 }
  element_counts.transform_keys! { |k| periodic_table_elements[k] } # vhange atomic mass to element symbol
  element_counts.transform_values! { |v| v > 1 ? v : '' } # remove 1s since, counts of 1 are implicit
end

def write_sum_formula_string(molecule, periodic_table_elements)
  # Write sum formula in the order C > H > all other elements in alphabetic order.
  element_counts = compute_element_counts(molecule, periodic_table_elements)
  sum_formula_string = ''
  sum_formula_string += "C#{element_counts['H']}" if element_counts.key?('C')
  sum_formula_string += "H#{element_counts['H']}" if element_counts.key?('H')
  element_counts.each do |element, count|
    sum_formula_string += "#{element}#{count}" unless %w[C H].include?(element)
  end
  puts "The sum formula of this molecule is: #{sum_formula_string}\n"
end

def write_ninchi_string(molecule, periodic_table_elements)
  ninchi_string = "nInChI=1S/#{write_sum_formula_string(molecule, periodic_table_elements)}/c#{serialization(molecule)}"
  puts "---------------------------------------"
  puts ninchi_string
  puts "---------------------------------------"
  ninchi_string
end

def canonicalize_molecule(molecule)
  fail "Structure is empty\n" if molecule.empty?

  puts "Initial molecule:"
  molecule.each { |atom| puts atom.inspect }

  sorted_molecule = sort_across_elements_by_atomic_mass(molecule)
  puts "Molecule with elements sorted by atomic mass (increasing):"
  sorted_molecule.each { |atom| puts atom.inspect }

  sorted_molecule = update_molecule_indices(sorted_molecule)
  puts "Molecule with updated indices after sorting:"
  sorted_molecule.each { |atom| puts atom.inspect }

  sorted_molecule = sort_connections_by_atomic_mass(sorted_molecule)
  puts "Molecule with connections sorted by atomic mass (decreasing):"
  sorted_molecule.each { |atom| puts atom.inspect }

  swap_logic1(sorted_molecule)
  puts "Molecule with connections sorted by atomic mass (decreasing):"
  sorted_molecule.each { |atom| puts atom.inspect }

  sorted_molecule = update_molecule_indices(sorted_molecule)
  puts "Molecule with updated indices after sorting:"
  sorted_molecule.each { |atom| puts atom.inspect }

  swap_logic2(sorted_molecule)
  puts "Molecule with connections sorted by atomic mass (decreasing):"
  sorted_molecule.each { |atom| puts atom.inspect }

  sorted_molecule = update_molecule_indices(sorted_molecule)
  puts "Molecule with updated indices after sorting:"
  sorted_molecule.each { |atom| puts atom.inspect }

  sorted_molecule
end

def serialization(molecule)
  raise '\nStructure is empty\n' if molecule.empty?

  graph = compute_graph(molecule)

  inchi_string = ''
  graph.each do |line|
    inchi_string += "(#{line[0]}-#{line[1]})" if line[0] != line[1]
  end
  inchi_string
end

def create_dot_file(molecule, periodic_table_colors)
  raise 'Structure is empty.' if molecule.empty?

  graph = compute_graph(molecule)

  dotfile = ''
  dotfile += "graph test\n"
  dotfile += "{\n"
  dotfile += "  bgcolor=grey\n"
  molecule.each_with_index do |atom, i|
    color = periodic_table_colors.fetch(atom[0], 'lightgrey')
    dotfile += "  #{i} [label=\"#{atom[0]} #{i}\" color=#{color},style=filled,shape=circle,fontname=Calibri];\n"
  end
  graph.each do |line|
    dotfile += "  #{line[0]} -- #{line[1]} [color=black,style=bold];\n" if line[0] != line[1]
  end
  dotfile += "}\n"

  puts "\nPrinting Connection Table\n"
  puts "Now printing graph\n"
  puts "---------------------------------------"
  puts dotfile
  puts "---------------------------------------"

  dotfile
end

# Helper methods ###################################################################################
private

def sort_connections_by_atomic_mass(molecule)
  # Sort connection numbers of each atom left to right from large to small.
  # Mutates `molecule`.
  sorted_molecule = []
  molecule.each do |atom|
    element, connections = atom
    sorted_connections = connections.sort { |a, b| b[1] <=> a[1] }
    sorted_molecule.push([element, sorted_connections])
  end
  sorted_molecule
end

def swap_logic1(molecule)
  # Mutates `molecule`.
  molecule.sort! do |atom_a, atom_b|
    case
    when (atom_a[0][1] == atom_b[0][1]) && # A and B are the same element ...
         (atom_a[1].length == atom_b[1].length) && # with the same number of connections ...
         (atom_a[1] <=> atom_b[1]) == 1 # and the connection indices of A are larger than the ones of B.
      1
    when (atom_a[0][1] == atom_b[0][1]) && # A and B are the same element ...
         (atom_a[1].length == atom_b[1].length) && # with the same number of connections ...
         (atom_a[1] <=> atom_b[1]) == -1 # and the connection indices of A are smaller than the ones of B.
      -1
    else
      0
    end
  end
end

def swap_logic2(molecule)
  # Mutates `molecule`.
  molecule.sort! do |atom_a, atom_b|
    case
    when (atom_a[0][1] == atom_b[0][1]) && # A and B are the same element ...
         (atom_a[1].length > atom_b[1].length) # and A has more connections than B.
      1
    when (atom_a[0][1] == atom_b[0][1]) && # A and B are the same element ...
         (atom_a[1].length < atom_b[1].length) # and A has less connections than B.
      -1
    else
      0
    end
  end
end

def sort_across_elements_by_atomic_mass(molecule)
  molecule.sort { |a, b| a[0][1] <=> b[0][1] }
end

def update_molecule_indices(molecule)
  index_updates = compute_index_updates(molecule)
  updated_molecule = []
  molecule.each do |atom|
    element, connections = atom
    updated_element = update_element_index(element, index_updates)
    updated_connections = update_connection_indices(connections, index_updates)
    updated_molecule.push([updated_element, updated_connections])
  end
  updated_molecule
end

def update_element_index(element, index_updates)
  element[0] = index_updates[element[0]] if index_updates.key?(element[0])
  element
end

def update_connection_indices(connections, index_updates)
  connections.map do |connection|
    connection[0] = index_updates[connection[0]] if index_updates.key?(connection[0])
    connection
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

def compute_graph(molecule)
  temp_molecule = Marshal.load(Marshal.dump(molecule))
  graph = []
  temp_molecule.each_with_index do |atom, i|
    atom.shift
    atom.each do |connected_atom|
      if i < connected_atom
        graph.push([i, connected_atom])
      else
        graph.push([connected_atom, i])
      end
    end
  end
  graph = graph.uniq.sort!
end

# def canonicalize_molecule(molecule, periodic_table_elements)
#   return molecule unless molecule.length > 1

#   previous_molecule_states = Set[molecule]
#   (0..(molecule.length - 1) * 20).each do |i|
#     puts "=== Start Pass ##{i} ===\n"
#     molecule = canonicalization(molecule, method(:swap_logic1), periodic_table_elements)
#     molecule = canonicalization(molecule, method(:swap_logic2), periodic_table_elements)
#     puts "=== End Pass ##{i} ===\n"
#     # Stop iterating if this molecule state occurred before. Recurrence of molecule states
#     # heuristically implies convergence in the form of either a) cyclic recurrence of
#     # molecule states (e.g., A, B, A) or b) unchanging/stable molecule state (e.g., A, A).
#     break if previous_molecule_states.include?(molecule)
#     previous_molecule_states.add(molecule)
#   end
#   molecule
# end

# def canonicalization(old_molecule, swap_logic, periodic_table_elements)
#   fail "Structure is empty\n" if old_molecule.empty?

#   puts "Now sorting the array\n"
#   puts "Old array:\n#{old_molecule}\n"

#   correspondence_table = compute_correspondence_table(periodic_table_elements, old_molecule)
#   new_molecule = sort_across_elements_by_atomic_mass(correspondence_table, old_molecule)
#   update_element_connection_indices(correspondence_table, new_molecule)
#   puts "New array withOUT labels re-organized:\n#{new_molecule}\n"

#   sort_within_elements_by_connection_indices(new_molecule)

#   atom_update = compute_atom_swaps(new_molecule, swap_logic)
#   return new_molecule unless atom_update

#   old_atom, new_atom = atom_update
#   new_molecule[old_atom], new_molecule[new_atom] = new_molecule[new_atom], new_molecule[old_atom]
#   swap_atoms(new_molecule, old_atom, new_atom)
#   sort_within_elements_by_connection_indices(new_molecule)
#   puts "\nSuggestion to re-order: #{new_molecule}\n"
#   new_molecule
#   # puts "New array WITH labels re-organized:\n#{new_molecule}\n"
# end

# def swap_logic1(molecule, index)
#   # Swap element A (molecule[index]) and element B (molecule[index + 1]) if ...
#   if (molecule[index][0] == molecule[index + 1][0]) && # A and B are the same element ...
#      (molecule[index].length == molecule[index + 1].length) && # with the same number of connections ...
#      (molecule[index] <=> molecule[index + 1]) == 1 # and the connection indices of A are larger than the ones of B.
#     old_atom = index
#     new_atom = index + 1
#     return [old_atom, new_atom]
#   end
#   nil
# end

# def swap_logic2(molecule, index)
#   # Swap element A (molecule[index]) and element B (molecule[index + 1]) if ...
#   if (molecule[index][0] == molecule[index + 1][0]) && # A and B are the same element ...
#      (molecule[index].length > molecule[index + 1].length) # and A has more connections than B.
#     old_atom = index
#     new_atom = index + 1
#     return [old_atom, new_atom]
#   end
#   nil
# end

# def compute_correspondence_table(periodic_table_elements, molecule)
#   # Compute correspondence hash table mapping old to new array positions.
#   correspondence_table = {}
#   molecule.each_with_index do |atom, i|
#     periodic_table_index = periodic_table_elements.index(atom[0])
#     correspondence_table[i] = periodic_table_index # map old index to index in periodic table
#   end
#   correspondence_table = correspondence_table.sort_by(&:last).to_h # sort by index in periodic table
#   correspondence_table.transform_values!.with_index { |_, idx| idx } # transform index in periodic table to new index
# end

# def compute_atom_swaps(molecule, swap_logic)
#   # If there's a swap while `atom_count` has not been reached, break and return
#   # `old_atom` and `new_atom`. Otherwise, if no swap occurs while iterating over
#   # atoms in `molecule` return nil.
#   last_atom_index = molecule.length - 1
#   puts "last_atom_index: #{last_atom_index} \n"
#   (0..last_atom_index - 1).each do |i|
#     puts "#{i}: #{molecule[i]}, #{i + 1}: #{molecule[i + 1]}"
#     swap_update = swap_logic.call(molecule, i)
#     next unless swap_update

#     old_atom, new_atom = swap_update
#     puts "Swap #{old_atom} vs. #{new_atom} \n"
#     return [old_atom, new_atom]
#   end
#   nil
# end

# def swap_atoms(molecule, old_atom, new_atom)
#   # Mutates `molecule`.
#   molecule.each do |element|
#     print "#{element} "
#     (1..element.length - 1).each do |i|
#       print "#{i}:#{element[i]}"
#       if element[i] == old_atom
#         print "c#{new_atom}"
#         element[i] = new_atom
#       elsif element[i] == new_atom
#         print "c#{old_atom}"
#         element[i] = old_atom
#       end
#       print ' '
#     end
#     print "\n"
#   end
# end