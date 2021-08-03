#
# (c) CC BY-SA | Ulrich Schatzschneider | Universität Würzburg | NFDI4Chem | v1.4 | 06.06.2021
#

require './periodic_table'
require 'set'

def read_molfile(filename)
  filename = 'NONE' if filename.empty?
  print "\nSupplied argument was: ", filename, "\n\n"
  if !File.exist?(filename)
    print 'File ', filename, " does not exist\n\n"
    return nil
  else
    print 'File ', filename, " exists\n\n"
  end
  File.open(filename, 'rt').read # open text file in read-only mode
end

def create_molecule_array(file, filename)
  input = []
  file.each_line do |line|
    input.push(line)
  end
  print "Printing molfile\n\n"
  print 'Start of file: ', filename, "\n"
  print "-----------------------------------\n"
  print file
  print "\n-----------------------------------\n"
  print 'End of file: ', filename, "\n\n"
  #
  # header of a molfile is 4 lines long
  #
  print "Printing molfile header\n"
  print "-----------------------------------\n"
  (0..3).each do |i|
    printf('%i: %s', i, input[i])
  end
  print "-----------------------------------\n"
  #
  # number of atoms in file is 1st number of 4th line (with index "3" as counting starts at zero)
  # number of bonds is 2nd number of 4th line (NOT: length of file minus header length (4 lines) minus number of atom definition lines { minus one (as final line is "M END") } )
  # (This is not true, since there can be additional lines with definitions starting with letter "M", thus rather use second number which is the number of bonds)
  #
  atom_count, bondCount = input[3].scan(/\d+/).map { |n| n.to_i }
  print "\nLength of file is ", input.length, " lines\n\n"
  print 'Number of atoms is ', atom_count, "\n\n"
  print 'Number of bonds is ', bondCount, "\n\n"
  #
  # now read the next lines containing the atom definitions
  # first three "fields" are "pseudo-coordinates", the 4th (with index 3 as counting starts at zero) is the element symbol which is what we want here, everything else is ignored
  #
  molecule = []
  (4..atom_count + 3).each do |i|
    atom = input[i].split(' ')
    printf("Line %i Atom %i: %s %s\n", i, i - 4, atom[3], atom)
    molecule.push([atom[3]])
  end
  print "\n"
  #
  # now read the remaining lines containing the bond definitions in the sequence atom1 atom2 bond_order ... unknown/unused ... (can be ignored)
  #
  (0..bondCount - 1).each do |i|
    connectionTable = input[i + 4 + atom_count].split(' ')
    if connectionTable[0] > connectionTable[1]
      connectionTable[0], connectionTable[1] = connectionTable[1], connectionTable[0] # make sure first atom always has lower (not: higher?) index
    end
    connectionTable[0] = String(connectionTable[0].to_i - 1)
    connectionTable[1] = String(connectionTable[1].to_i - 1)
    printf("Bond %i: %s\n", i + 1, connectionTable[0] + '-' + connectionTable[1])
    molecule[connectionTable[0].to_i].push(connectionTable[1].to_i)    # need to push twice, to the first atom of a bond
    molecule[connectionTable[1].to_i].push(connectionTable[0].to_i)    # and then to the second atom of the bond
  end
  sortedMolecule = []
  molecule.each do |atom|
    atomSorted = []
    elementSymbol = atom[0]
    atom.shift
    atomSorted = atom.sort.reverse
    atom = atomSorted
    atom.unshift(elementSymbol)
    sortedMolecule.push(atom)
  end
  molecule = sortedMolecule
  #
  # The "real" routine ends here, the rest is just around for test purpose
  #
  print "\nPrinting Connection Table\n\n"
  #
  # print the "connection table" in format "atom number: element_symbol connected_atoms"
  # start with highest priority atom
  #
  molecule.reverse_each do |atom|
    print atom_count - 1, ':'
    atom.each do |connectionTable|
      print ' ' + connectionTable.to_s
    end
    atom_count -= 1
    print "\n"
  end
  print "\n", molecule, "\n\n"
  molecule
end

def calculate_sum_formula(molecule)
  periodic_table_elements = []
  periodic_table_elements = PeriodicTable::Elements
  #
  # create sum formula array in the order C > H > all other elements in alphabetic order and "atom count" set to zero
  #
  sumFormula = []
  sumFormula.push(['C', 0], ['H', 0])
  periodic_table_elements = periodic_table_elements - ['C'] - ['H']
  periodic_table_elements.sort!
  periodic_table_elements.each do |element|
    sumFormula.push([element, 0])
  end
  #
  # sum up
  #
  sumFormula.each do |element|
    molecule.each do |atom|
      element[1] = element[1] + 1 if element[0] == atom[0]
    end
  end
  #
  # generate output string
  #
  sumFormulaString = ''
  sumFormula.each do |element|
    if element[1] > 1
      sumFormulaString = sumFormulaString + element[0].to_s + String(element[1])
    elsif element[1] > 0
      sumFormulaString += element[0].to_s # if there is only one atom of a particular element in the molecule, only add the element symbol without stoichiometric count
    end
  end
  sumFormulaString = 'ERROR - This *.mol file does not contain any molecular structure' if sumFormulaString == ''
  #
  # print sum formula
  #
  print 'The sum formula of this molecule is: ', sumFormulaString, "\n\n"
  print sumFormula, "\n\n"
  sumFormulaString
end

def serialization(molecule)
  if molecule == []
    print "\nStructure is empty\n"
    return 'Structure is empty'
  else
    print molecule, "\n\n"
  end
  tempMolecule = []
  tempMolecule = Marshal.load(Marshal.dump(molecule))
  graph = []
  i = 0
  tempMolecule.each do |atom|
    atom.shift
    atom.each do |connectedAtom|
      if i < connectedAtom
        graph.push([i, connectedAtom])
      else
        graph.push([connectedAtom, i])
      end
    end
    i += 1
  end
  graph = graph.uniq.sort!
  print "Now printing new InChI\n\n"
  print "---------------------------------------\n"
  print 'InChI=1S/sum_formula/'
  graph.each do |line|
    print '(', line[0], '-', line[1], ')' if line[0] != line[1]
  end
  print "\n---------------------------------------\n"
  #
  # add output to inChIstring
  #
  inChIstring = ''
  graph.each do |line|
    inChIstring = inChIstring + '(' + line[0].to_s + '-' + line[1].to_s + ')' if line[0] != line[1]
  end
  inChIstring
end

def create_dot_file(molecule)
  dotfile = ''
  print "\nPrinting Connection Table\n\n"
  if molecule == []
    print "Structure is empty\n"
    return 'Structure is empty'
  else
    print molecule, "\n\n"
  end
  tempMolecule = []
  tempMolecule = Marshal.load(Marshal.dump(molecule))
  graph = []
  i = 0
  tempMolecule.each do |atom|
    atom.shift
    atom.each do |connectedAtom|
      if i < connectedAtom
        graph.push([i, connectedAtom])
      else
        graph.push([connectedAtom, i])
      end
    end
    i += 1
  end
  graph = graph.uniq.sort!
  elementColor = PeriodicTable::ElementColor
  #
  # print output to console
  #
  print "Now printing graph\n\n"
  print "---------------------------------------\n"
  print "graph test\n"
  print "{\n"
  print "  bgcolor=grey\n"
  i = 0
  molecule.each do |atom|
    color = elementColor.fetch(atom[0], 'lightgrey') # if color is unspecified fall back on lightgray
    print '  ', i, ' [label="', atom[0], ' ', i, '" color=', color, ",style=filled,shape=circle,fontname=Calibri];\n"
    i += 1
  end
  graph.each do |line|
    print '  ', line[0], ' -- ', line[1], " [color=black,style=bold];\n" if line[0] != line[1]
  end
  print "}\n"
  print "---------------------------------------\n"
  #
  # add output to dotfile
  #
  dotfile += "graph test\n"
  dotfile += "{\n"
  dotfile += "  bgcolor=grey\n"
  i = 0
  molecule.each do |atom|
    color = elementColor.fetch(atom[0], 'lightgrey')
    dotfile = dotfile + '  ' + i.to_s + ' [label="' + atom[0].to_s + ' ' + i.to_s + '" color=' + color + ",style=filled,shape=circle,fontname=Calibri];\n"
    i += 1
  end
  graph.each do |line|
    if line[0] != line[1]
      dotfile = dotfile + '  ' + line[0].to_s + ' -- ' + line[1].to_s + " [color=black,style=bold];\n"
    end
  end
  dotfile + "}\n"
end

def create_ninchi_string(molecule)
  "nInChI=1S/#{calculate_sum_formula(molecule)}/c#{serialization(molecule)}"
end

def sort_connection_numbers(molecule)
  # Sort connection numbers of each atom left to right from large to small.
  # Mutates `molecule`.
  molecule.each do |atom|
    element_symbol = atom[0]
    atom.shift
    atom.sort!.reverse! if atom.length > 1
    atom.insert(0, element_symbol)
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

def canonicalization(old_molecule, swap_logic, periodic_table_elements)
  fail "Structure is empty\n" if old_molecule.empty?

  atom_count = old_molecule.length - 1 # determine number of rows, make sure to properly account for array starting at index zero

  puts "Now sorting the array\n"
  sort_connection_numbers(old_molecule)
  puts "Old array:\n#{old_molecule}\n"

  new_molecule = compute_new_molecule(periodic_table_elements, old_molecule, atom_count)
  puts "New array withOUT labels re-organized:\n#{new_molecule}\n"

  correspondence_table = compute_correspondance_table(periodic_table_elements, old_molecule, atom_count)

  old_molecule = Marshal.load(Marshal.dump(new_molecule)) # marshalling results in deep (as opposed to shallow) copy

  new_molecule = compute_element_connections(old_molecule, correspondence_table, atom_count)
  sort_connection_numbers(new_molecule)

  puts "atom_count: #{atom_count} \n"
  atom_update = compute_atom_swaps(new_molecule, atom_count, swap_logic)
  return new_molecule unless atom_update
  old_atom, new_atom = atom_update
  old_molecule = Marshal.load(Marshal.dump(new_molecule))
  old_molecule[old_atom], old_molecule[new_atom] = old_molecule[new_atom], old_molecule[old_atom]
  swap_atoms(old_molecule, old_atom, new_atom, atom_count)
  sort_connection_numbers(old_molecule)
  puts "\nSuggestion to re-order: #{old_molecule}\n"
  old_molecule
  # puts "New array WITH labels re-organized:\n#{new_molecule}\n"
end

def canonicalize_molecule(molecule)
  return molecule unless molecule.length > 1

  previous_molecule_states = Set[molecule]
  (0..(molecule.length - 1) * 20).each do |i|
    puts "=== Start Pass ##{i} ===\n"
    molecule = canonicalization(molecule, method(:swap_logic1), PeriodicTable::Elements)
    molecule = canonicalization(molecule, method(:swap_logic2), PeriodicTable::Elements)
    puts "=== End Pass ##{i} ===\n"
    # Stop iterating if this molecule state occurred before. Recurrence of molecule states
    # heuristically implies convergence in the form of either a) cyclic recurrence of
    # molecule states (e.g., A, B, A) or b) unchanging/stable molecule state (e.g., A, A).
    break if previous_molecule_states.include?(molecule)
    previous_molecule_states.add(molecule)
  end
  molecule
end

def compute_new_molecule(periodic_table_elements, molecule, atom_count)
  new_molecule = []
  periodic_table_elements.each do |element| # second pass - sort by element in increasing order, lowest atomic mass element to the left/bottom
    (0..atom_count).each do |i|
      next unless molecule[i][0] == element

      new_molecule.push(molecule[i])
    end
  end
  new_molecule
end

def compute_correspondance_table(periodic_table_elements, molecule, atom_count)
  # Compute correspondance table containing pairs of array positions [old,new].
  correspondence_table = []
  j = 0
  periodic_table_elements.each do |element| # sort by element in increasing order, lowest atomic mass element to the left/bottom
    (0..atom_count).each do |i|
      next unless molecule[i][0] == element

      correspondence_table.push([i, j])
      j += 1
    end
  end
  correspondence_table.sort! { |a, b| b[0] <=> a[0] } # sort (in-place) by lowest old atom position
end

def compute_element_connections(molecule, correspondence_table, atom_count)
  new_molecule = []
  (0..atom_count).each do |i|
    temp_array = []
    temp_array.push(molecule[i][0]) # add element symbol to temporary array
    (1..molecule[i].length - 1).each do |j|
      (0..correspondence_table.length - 1).each do |k|
        if molecule[i][j] == correspondence_table[k][0]
          temp_array.push(correspondence_table[k][1]) # append new connection to temporary array
        end
      end
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
