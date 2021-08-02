#
# (c) CC BY-SA | Ulrich Schatzschneider | Universität Würzburg | NFDI4Chem | v1.4 | 06.06.2021
#

require './periodic_table'

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
  atomCount, bondCount = input[3].scan(/\d+/).map { |n| n.to_i }
  print "\nLength of file is ", input.length, " lines\n\n"
  print 'Number of atoms is ', atomCount, "\n\n"
  print 'Number of bonds is ', bondCount, "\n\n"
  #
  # now read the next lines containing the atom definitions
  # first three "fields" are "pseudo-coordinates", the 4th (with index 3 as counting starts at zero) is the element symbol which is what we want here, everything else is ignored
  #
  molecule = []
  (4..atomCount + 3).each do |i|
    atom = input[i].split(' ')
    printf("Line %i Atom %i: %s %s\n", i, i - 4, atom[3], atom)
    molecule.push([atom[3]])
  end
  print "\n"
  #
  # now read the remaining lines containing the bond definitions in the sequence atom1 atom2 bond_order ... unknown/unused ... (can be ignored)
  #
  (0..bondCount - 1).each do |i|
    connectionTable = input[i + 4 + atomCount].split(' ')
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
    print atomCount - 1, ':'
    atom.each do |connectionTable|
      print ' ' + connectionTable.to_s
    end
    atomCount -= 1
    print "\n"
  end
  print "\n", molecule, "\n\n"
  molecule
end

def calculate_sum_formula(molecule)
  periodicTable = []
  periodicTable = PeriodicTable::Elements
  #
  # create sum formula array in the order C > H > all other elements in alphabetic order and "atom count" set to zero
  #
  sumFormula = []
  sumFormula.push(['C', 0], ['H', 0])
  periodicTable = periodicTable - ['C'] - ['H']
  periodicTable.sort!
  periodicTable.each do |element|
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

def canonicalization1(old_molecule)
  if old_molecule == []
    print "Structure is empty\n"
    return 'Structure is empty'
  end
  periodicTable = PeriodicTable::Elements
  atomCount = old_molecule.length - 1 # determine number of rows, make sure to properly account for array starting at index zero
  printf("Now sorting the array\n\n")
  sort_connection_numbers(old_molecule)
  print "Old array: \n\n", old_molecule, "\n\n"
  new_molecule = []
  correspondenceTable = [] # contains pairs of array positions [old,new]
  j = 0
  periodicTable.each do |element| # second pass - sort by element in increasing order, lowest atomic mass element to the left/bottom
    (0..atomCount).each do |i|
      next unless old_molecule[i][0] == element

      new_molecule.push(old_molecule[i])
      correspondenceTable.push([i, j])
      j += 1
    end
  end
  old_molecule = Marshal.load(Marshal.dump(new_molecule))
  temp = # sort the correspondence table by lowest old atom position
    correspondenceTable.sort do |a, b|
      b[0] <=> a[0]
    end
  correspondenceTable = temp
  print "New array withOUT labels re-organized: \n\n", new_molecule, "\n\n"
  new_molecule = []
  (0..atomCount).each do |i|
    tempArray = []
    tempArray.push(old_molecule[i][0]) # add element symbol to temporary array
    (1..old_molecule[i].length - 1).each do |j|
      (0..correspondenceTable.length - 1).each do |k|
        if old_molecule[i][j] == correspondenceTable[k][0]
          tempArray.push(correspondenceTable[k][1]) # append new connection to temporary array
        end
      end
    end
    new_molecule.push(tempArray) # add whole new atom connection list for atom to temporary array
  end
  sort_connection_numbers(new_molecule)
  print 'atomCount: ', atomCount, "\n\n"
  i = 0
  exitLoop = false
  until exitLoop
    print i, ': ', new_molecule[i], ' ', i + 1, ': ', new_molecule[i + 1], "\n"
    if new_molecule[i][0] == new_molecule[i + 1][0] && (new_molecule[i].length > new_molecule[i + 1].length)
      oldAtom = i
      newAtom = i + 1
      print 'Swap ', oldAtom, ' vs. ', newAtom, "\n\n"
      exitLoop = true
    end
    i += 1
    next unless i >= atomCount

    exitLoop = true
    print "\n"
    return new_molecule
  end
  old_molecule = Marshal.load(Marshal.dump(new_molecule))
  old_molecule[oldAtom], old_molecule[newAtom] = old_molecule[newAtom], old_molecule[oldAtom]
  (0..atomCount).each do |i|
    line = []
    line = old_molecule[i]
    print line, ' '
    (1..line.length - 1).each do |j|
      print j, ':', line[j]
      if line[j] == oldAtom
        print 'c', newAtom
        line[j] = newAtom
      elsif line[j] == newAtom
        print 'c', oldAtom
        line[j] = oldAtom
      end
      print ' '
    end
    print "\n"
  end
  sort_connection_numbers(old_molecule)
  print "\nSuggestion to re-order: ", old_molecule, "\n\n"
  old_molecule
  #
  #
  #
  #    print "CorrespondenceTable [old,new]: \n\n",correspondenceTable,"\n\n"
  #    print "New array WITH labels re-organized: \n\n",new_molecule,"\n\n"
  #    return new_molecule
end

def canonicalization2(old_molecule)
  if old_molecule == []
    print "Structure is empty\n"
    return 'Structure is empty'
  end
  periodicTable = PeriodicTable::Elements
  atomCount = old_molecule.length - 1 # determine number of rows, make sure to properly account for array starting at index zero
  printf("Now sorting the array\n\n")
  sort_connection_numbers(old_molecule)
  print "Old array: \n\n", old_molecule, "\n\n"
  new_molecule = []
  correspondenceTable = [] # contains pairs of array positions [old,new]
  j = 0
  periodicTable.each do |element| # second pass - sort by element in increasing order, lowest atomic mass element to the left/bottom
    (0..atomCount).each do |i|
      next unless old_molecule[i][0] == element

      new_molecule.push(old_molecule[i])
      correspondenceTable.push([i, j])
      j += 1
    end
  end
  old_molecule = Marshal.load(Marshal.dump(new_molecule))
  temp = # sort the correspondence table by lowest old atom position
    correspondenceTable.sort do |a, b|
      b[0] <=> a[0]
    end
  correspondenceTable = temp
  print "New array withOUT labels re-organized: \n\n", new_molecule, "\n\n"
  new_molecule = []
  (0..atomCount).each do |i|
    tempArray = []
    tempArray.push(old_molecule[i][0]) # add element symbol to temporary array
    (1..old_molecule[i].length - 1).each do |j|
      (0..correspondenceTable.length - 1).each do |k|
        if old_molecule[i][j] == correspondenceTable[k][0]
          tempArray.push(correspondenceTable[k][1]) # append new connection to temporary array
        end
      end
    end
    new_molecule.push(tempArray) # add whole new atom connection list for atom to temporary array
  end
  sort_connection_numbers(new_molecule)
  print 'atomCount: ', atomCount, "\n\n"
  i = 0
  exitLoop = false
  until exitLoop
    print i, ': ', new_molecule[i], ' ', i + 1, ': ', new_molecule[i + 1], "\n"
    if (new_molecule[i][0] == new_molecule[i + 1][0]) && (new_molecule[i].length == new_molecule[i + 1].length)
      do_swap = (new_molecule[i] <=> new_molecule[i + 1])
      #        do_swap=(new_molecule[i].to_s <=> new_molecule[i+1].to_s)
      #        if(new_molecule[i].to_s > new_molecule[i+1].to_s)
      if do_swap == 1
        oldAtom = i
        newAtom = i + 1
        print 'Swap ', oldAtom, ' vs. ', newAtom, "\n\n"
        exitLoop = true
      end
    end
    i += 1
    next unless i >= atomCount

    exitLoop = true
    print "\n"
    return new_molecule
  end
  old_molecule = Marshal.load(Marshal.dump(new_molecule))
  old_molecule[oldAtom], old_molecule[newAtom] = old_molecule[newAtom], old_molecule[oldAtom]
  (0..atomCount).each do |i|
    line = []
    line = old_molecule[i]
    print line, ' '
    (1..line.length - 1).each do |j|
      print j, ':', line[j]
      if line[j] == oldAtom
        print 'c', newAtom
        line[j] = newAtom
      elsif line[j] == newAtom
        print 'c', oldAtom
        line[j] = oldAtom
      end
      print ' '
    end
    print "\n"
  end
  sort_connection_numbers(old_molecule)
  print "\nSuggestion to re-order: ", old_molecule, "\n\n"
  old_molecule
  #
  #
  #
  #    print "CorrespondenceTable [old,new]: \n\n",correspondenceTable,"\n\n"
  #    print "New array WITH labels re-organized: \n\n",new_molecule,"\n\n"
  #    return new_molecule
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

def canonicalize_molecule(molecule)
  if molecule.length > 1
    (0..(molecule.length - 1) * 20).each do |i|
      print '=== Start Pass #', i, " ===\n\n"
      molecule = canonicalization1(molecule)
      molecule = canonicalization2(molecule)
      print '=== End Pass #', i, " ===\n\n"
    end
  end
  molecule
end

def create_ninchi_string(molecule)
  sumFormulaString = calculate_sum_formula(molecule)
  "nInChI=1S/#{sumFormulaString}/c#{serialization(molecule)}"
end

def sort_connection_numbers(molecule)
  # Sort connection numbers of each atom left to right
  # from large to small.
  # Mutates `molecule`.
  molecule.each do |atom|
    element_symbol = atom[0]
    atom.shift
    atom.sort!.reverse! if atom.length > 1
    atom.insert(0, element_symbol)
  end
end
