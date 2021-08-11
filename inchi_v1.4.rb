# (c) CC BY-SA | Ulrich Schatzschneider | Universität Würzburg | NFDI4Chem | v1.4 | 06.06.2021

require 'set'

def read_molfile(filename)
  raise "#{filename} doesn't exist." unless File.exist?(filename)

  molfile = File.read(filename) # reads entire file and closes it
  molfile_lines = molfile.split("\n")

  puts "Printing molfile\n"
  puts "Start of file: #{filename}\n"
  puts "-----------------------------------"
  puts molfile
  puts "-----------------------------------"
  puts "End of file: #{filename}\n"

  puts "Printing molfile header"
  puts "-----------------------------------"
  (0..3).each do |i| # header of a molfile is 4 lines long
    puts "#{i}: #{molfile_lines[i]}\n"
  end
  puts "-----------------------------------"

  molfile_lines
end

def create_molecule_array(molfile_lines)
  # number of atoms in file is 1st number of 4th line (with index "3" as counting starts at zero)
  # number of bonds is 2nd number of 4th line (NOT: length of file minus header length (4 lines) minus number of atom definition lines { minus one (as final line is "M END") } )
  # (This is not true, since there can be additional lines with definitions starting with letter "M", thus rather use second number which is the number of bonds)
  atom_count, bond_count = molfile_lines[3].scan(/\d+/).map { |n| n.to_i }
  puts "\nLength of file is #{molfile_lines.length} lines\n"
  puts "Number of atoms is #{atom_count}\n"
  puts "Number of bonds is #{bond_count}\n"

  # now read the next lines containing the atom definitions
  # first three "fields" are "pseudo-coordinates", the 4th (with index 3 as counting starts at zero) is the element symbol which is what we want here, everything else is ignored
  molecule = []
  (4..atom_count + 3).each do |atom_index|
    atom = molfile_lines[atom_index].split(' ')
    puts "Line #{atom_index} Atom #{atom_index - 4}: #{atom[3]} #{atom}"
    molecule.push([atom[3]])
  end

  # now read the remaining lines containing the bond definitions in the sequence atom1 atom2 bond_order ... unknown/unused ... (can be ignored)
  (0..bond_count - 1).each do |bond_index|
    vertex1, vertex2, * = molfile_lines[bond_index + 4 + atom_count].split(' ').map { |i| i.to_i - 1 }
    vertex1, vertex2 = vertex2, vertex1 if vertex1 > vertex2 # make sure first atom always has lower (not: higher?) index
    puts "Bond #{bond_index + 1}: #{vertex1}-#{vertex2}"
    molecule[vertex1].push(vertex2)    # need to push twice, to the first atom of a bond
    molecule[vertex2].push(vertex1)    # and then to the second atom of the bond
  end

  sort_connection_numbers(molecule) # method ends here
  # puts molecule.inspect # printing connection table is just around for test purpose
end

def calculate_sum_formula(molecule, periodic_table_elements)
  # Create sum formula array in the order C > H > all other elements in alphabetic order.
  # "atom count" is set to zero.
  sum_formula = []
  sum_formula.push(['C', 0], ['H', 0])
  periodic_table_elements = periodic_table_elements - ['C'] - ['H']
  periodic_table_elements.sort!
  periodic_table_elements.each do |element|
    sum_formula.push([element, 0])
  end
  sum_formula.each do |element|
    molecule.each do |atom|
      element[1] += 1 if element[0] == atom[0] # sum up
    end
  end

  # generate output string
  sum_formula_string = ''
  sum_formula.each do |element|
    if element[1] > 1
      sum_formula_string = sum_formula_string + element[0].to_s + String(element[1])
    elsif element[1] > 0
      sum_formula_string += element[0].to_s # if there is only one atom of a particular element in the molecule, only add the element symbol without stoichiometric count
    end
  end
  raise 'ERROR - This *.mol file does not contain any molecular structure' if sum_formula_string == ''

  puts "The sum formula of this molecule is: #{sum_formula_string}\n"
  print sum_formula, "\n\n"
  sum_formula_string
end

def canonicalization(old_molecule, swap_logic, periodic_table_elements)
  fail "Structure is empty\n" if old_molecule.empty?

  atom_count = old_molecule.length - 1 # determine number of rows, make sure to properly account for array starting at index zero

  puts "Now sorting the array\n"
  sort_connection_numbers(old_molecule)
  puts "Old array:\n#{old_molecule}\n"

  new_molecule = sort_elements_by_atomic_mass(periodic_table_elements, old_molecule)
  puts "New array withOUT labels re-organized:\n#{new_molecule}\n"

  correspondence_table = compute_correspondance_table(periodic_table_elements, old_molecule)

  new_molecule = compute_element_connections(new_molecule, correspondence_table)
  sort_connection_numbers(new_molecule)

  puts "atom_count: #{atom_count} \n"
  atom_update = compute_atom_swaps(new_molecule, atom_count, swap_logic)
  return new_molecule unless atom_update

  old_atom, new_atom = atom_update
  new_molecule[old_atom], new_molecule[new_atom] = new_molecule[new_atom], new_molecule[old_atom]
  swap_atoms(new_molecule, old_atom, new_atom, atom_count)
  sort_connection_numbers(new_molecule)
  puts "\nSuggestion to re-order: #{new_molecule}\n"
  new_molecule
  # puts "New array WITH labels re-organized:\n#{new_molecule}\n"
end

def serialization(molecule)
  raise '\nStructure is empty\n' if molecule.empty?

  graph = compute_graph(molecule)

  inchi_string = ''
  graph.each do |line|
    inchi_string += "(#{line[0]}-#{line[1]})" if line[0] != line[1]
  end
  puts "Now printing new InChI\n"
  puts "---------------------------------------"
  puts "InChI=1S/sum_formula/#{inchi_string}"
  puts "---------------------------------------"
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

def create_ninchi_string(molecule, periodic_table_elements)
  "nInChI=1S/#{calculate_sum_formula(molecule, periodic_table_elements)}/c#{serialization(molecule)}"
end

def canonicalize_molecule(molecule, periodic_table_elements)
  return molecule unless molecule.length > 1

  previous_molecule_states = Set[molecule]
  (0..(molecule.length - 1) * 20).each do |i|
    puts "=== Start Pass ##{i} ===\n"
    molecule = canonicalization(molecule, method(:swap_logic1), periodic_table_elements)
    molecule = canonicalization(molecule, method(:swap_logic2), periodic_table_elements)
    puts "=== End Pass ##{i} ===\n"
    # Stop iterating if this molecule state occurred before. Recurrence of molecule states
    # heuristically implies convergence in the form of either a) cyclic recurrence of
    # molecule states (e.g., A, B, A) or b) unchanging/stable molecule state (e.g., A, A).
    break if previous_molecule_states.include?(molecule)
    previous_molecule_states.add(molecule)
  end
  molecule
end

# Helper methods ###################################################################################S
private

def sort_connection_numbers(molecule)
  # Sort connection numbers of each atom left to right from large to small.
  # Mutates `molecule`.
  molecule.map! do |atom|
    element_symbol, *connections = atom
    connections.sort!.reverse! if connections.length > 1
    [element_symbol] + connections
  end
end

def swap_logic1(molecule, index)
  if (molecule[index][0] == molecule[index + 1][0]) && (molecule[index].length == molecule[index + 1].length)
    if (molecule[index] <=> molecule[index + 1]) == 1
      old_atom = index
      new_atom = index + 1
      return [old_atom, new_atom]
    end
  end
  nil
end

def swap_logic2(molecule, index)
  if molecule[index][0] == molecule[index + 1][0] && (molecule[index].length > molecule[index + 1].length)
    old_atom = index
    new_atom = index + 1
    return [old_atom, new_atom]
  end
  nil
end

def sort_elements_by_atomic_mass(periodic_table_elements, molecule)
  new_molecule = []
  periodic_table_elements.each do |element| # second pass - sort by element in increasing order, lowest atomic mass element to the left/bottom
    molecule.each do |molecule_element|
      next unless molecule_element[0] == element

      new_molecule.push(molecule_element)
    end
  end
  new_molecule
end

def compute_correspondance_table(periodic_table_elements, molecule)
  # Compute correspondance hash table mapping old to new array positions.
  correspondence_table = {}
  j = 0
  periodic_table_elements.each do |element| # sort by element in increasing order, lowest atomic mass element to the left/bottom
    molecule.each_with_index do |molecule_element, i|
      next unless molecule_element[0] == element

      correspondence_table[i] = j
      j += 1
    end
  end
  correspondence_table.sort.to_h # sort (in-place) by lowest old atom position
end

def compute_element_connections(molecule, correspondence_table)
  new_molecule = []
  molecule.each do |element|
    temp_array = []
    element_symbol, *connections = element
    temp_array.push(element_symbol)
    connections.each do |connection|
      temp_array.push(correspondence_table[connection])
    end
    new_molecule.push(temp_array) # add whole new atom connection list for atom to temporary array
  end
  new_molecule
end

def compute_atom_swaps(molecule, atom_count, swap_logic)
  # If there's a swap while `atom_count` has not been reached, break and return
  # `old_atom` and `new_atom`. Otherwise, if no swap occurs while iterating over
  # atoms in `molecule` return nil.
  (0..atom_count - 1).each do |i|
    puts "#{i}: #{molecule[i]}, #{i + 1}: #{molecule[i + 1]}"
    swap_update = swap_logic.call(molecule, i)
    next unless swap_update

    old_atom, new_atom = swap_update
    puts "Swap #{old_atom} vs. #{new_atom} \n"
    return [old_atom, new_atom]
  end
  nil
end

def swap_atoms(molecule, old_atom, new_atom, atom_count)
  # Mutates `molecule`.
  (0..atom_count).each do |i|
    print "#{molecule[i]} "
    (1..molecule[i].length - 1).each do |j|
      print "#{j}:#{molecule[i][j]}"
      if molecule[i][j] == old_atom
        print "c#{new_atom}"
        molecule[i][j] = new_atom
      elsif molecule[i][j] == new_atom
        print "c#{old_atom}"
        molecule[i][j] = old_atom
      end
      print ' '
    end
    print "\n"
  end
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
