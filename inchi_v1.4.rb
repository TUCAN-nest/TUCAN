#
# (c) CC BY-SA | Ulrich Schatzschneider | Universität Würzburg | NFDI4Chem | v1.4 | 06.06.2021
#

require './periodic_table'


def read_molfile(filename)
  filename='NONE' if filename.empty?
  print "\nSupplied argument was: ",filename,"\n\n"
  if(!File.exist?(filename))
    print "File ",filename," does not exist\n\n"
    return nil
  else
    print "File ",filename, " exists\n\n"
  end
  file=File.open(filename,"rt").read          # open text file in read-only mode
  return file
end


def create_molecule_array(file,filename)
  input=Array.new
  file.each_line do |line|
    input.push(line)
  end
  #
  print "Printing molfile\n\n"
  #
  print "Start of file: ",filename,"\n"
  print "-----------------------------------\n"
  print file
  print "\n-----------------------------------\n"
  print "End of file: ",filename,"\n\n"
  #
  # header of a molfile is 4 lines long
  #
  print "Printing molfile header\n"
  print "-----------------------------------\n"
  for i in 0..3
    printf("%i: %s",i,input[i])
  end
  print "-----------------------------------\n"
  #
  # number of atoms in file is 1st number of 4th line (with index "3" as counting starts at zero)
  # number of bonds is 2nd number of 4th line (NOT: length of file minus header length (4 lines) minus number of atom definition lines { minus one (as final line is "M END") } )
  # (This is not true, since there can be additional lines with definitions starting with letter "M", thus rather use second number which is the number of bonds)
  #
  atomCount, bondCount = input[3].scan(/\d+/).map { |n| n.to_i }
  #
  print "\nLength of file is ",input.length," lines\n\n"
  print "Number of atoms is ",atomCount,"\n\n"
  print "Number of bonds is ",bondCount,"\n\n"
  #
  # now read the next lines containing the atom definitions
  # first three "fields" are "pseudo-coordinates", the 4th (with index 3 as counting starts at zero) is the element symbol which is what we want here, everything else is ignored
  #
  molecule=Array.new
  for i in 4..atomCount+3
    atom=input[i].split(" ")
    printf("Line %i Atom %i: %s %s\n",i,i-4,atom[3],atom)
    molecule.push([atom[3]])
  end
  print "\n"
  #
  # now read the remaining lines containing the bond definitions in the sequence atom1 atom2 bond_order ... unknown/unused ... (can be ignored)
  #
  for i in 0..bondCount-1
    connectionTable=input[i+4+atomCount].split(" ")
    if connectionTable[0] > connectionTable[1]
      connectionTable[0],connectionTable[1] = connectionTable[1],connectionTable[0]  # make sure first atom always has lower (not: higher?) index
    end
    connectionTable[0]=String(connectionTable[0].to_i-1)
    connectionTable[1]=String(connectionTable[1].to_i-1)
    printf("Bond %i: %s\n",i+1,connectionTable[0]+"-"+connectionTable[1])
    molecule[connectionTable[0].to_i].push(connectionTable[1].to_i)    # need to push twice, to the first atom of a bond
    molecule[connectionTable[1].to_i].push(connectionTable[0].to_i)    # and then to the second atom of the bond
  end
  sortedMolecule=Array.new
  molecule.each do |atom|
    atomSorted=Array.new
    elementSymbol=atom[0]
    atom.shift
    atomSorted=atom.sort.reverse
    atom=atomSorted
    atom.unshift(elementSymbol)
    sortedMolecule.push(atom)
  end
  molecule=sortedMolecule
  #
  # The "real" routine ends here, the rest is just around for test purpose
  #
  print "\nPrinting Connection Table\n\n"
  #
  # print the "connection table" in format "atom number: element_symbol connected_atoms"
  # start with highest priority atom
  #
  molecule.reverse_each do |atom|
      print atomCount-1,":"
      atom.each do |connectionTable|
          print " "+connectionTable.to_s
      end
      atomCount=atomCount-1
      print "\n"
  end
  print "\n",molecule,"\n\n"
  return molecule
end


def calculate_sum_formula(molecule)
  periodicTable=Array.new
  periodicTable=PeriodicTable::Elements
  #
  # create sum formula array in the order C > H > all other elements in alphabetic order and "atom count" set to zero
  #
  sumFormula=Array.new
  sumFormula.push(['C',0],['H',0])
  periodicTable=periodicTable-['C']-['H']
  periodicTable.sort!
  periodicTable.each do |element|
    sumFormula.push([element,0])
  end
  #
  # sum up
  #
  sumFormula.each do |element|
      molecule.each do |atom|
          if(element[0] == atom[0])
              element[1]=element[1]+1
          end
      end
  end
  #
  # generate output string
  #
  sumFormulaString=''
  sumFormula.each do |element|
      if(element[1] > 1)
        sumFormulaString=sumFormulaString+element[0].to_s+String(element[1])
      elsif(element[1] > 0)
        sumFormulaString=sumFormulaString+element[0].to_s   # if there is only one atom of a particular element in the molecule, only add the element symbol without stoichiometric count
      end
  end
  if(sumFormulaString == '')
    sumFormulaString='ERROR - This *.mol file does not contain any molecular structure'
  end
  #
  # print sum formula
  #
  print "The sum formula of this molecule is: ",sumFormulaString,"\n\n"
  print sumFormula,"\n\n"
  return sumFormulaString
end


def canonicalization1(old_molecule)
  if(old_molecule == [])
    print "Structure is empty\n"
    return "Structure is empty"
  end
  periodicTable=Array.new
  periodicTable=PeriodicTable::Elements
  atomCount=old_molecule.length-1           # determine number of rows, make sure to properly account for array starting at index zero
  printf("Now sorting the array\n\n")
  old_molecule.each do |atom|   # first pass - sort connection numbers of each atom left to right from large to small
    elementSymbol=atom[0]
    atom.shift
    if(atom.length > 1)
      atom.sort!.reverse!
    end
    atom.insert(0,elementSymbol)
  end
  print "Old array: \n\n",old_molecule,"\n\n"
  new_molecule=Array.new
  correspondenceTable=Array.new   # contains pairs of array positions [old,new]
  j=0
  periodicTable.each do |element| # second pass - sort by element in increasing order, lowest atomic mass element to the left/bottom
    for i in 0..atomCount
      if(old_molecule[i][0] == element)
        new_molecule.push(old_molecule[i])
        correspondenceTable.push([i,j])
        j=j+1
      end
    end
  end
  old_molecule=Marshal.load(Marshal.dump(new_molecule))
  temp=Array.new
  temp=correspondenceTable.sort{|a,b| b[0] <=> a[0]}.reverse   # sort the correspondence table by lowest old atom position
  correspondenceTable=temp
  print "New array withOUT labels re-organized: \n\n",new_molecule,"\n\n"
  new_molecule=Array.new
  for i in 0..atomCount
    tempArray=Array.new
    tempArray.push(old_molecule[i][0])   # add element symbol to temporary array
    for j in 1..old_molecule[i].length-1
      for k in 0..correspondenceTable.length-1
        if(old_molecule[i][j] == correspondenceTable[k][0])
          tempArray.push(correspondenceTable[k][1])      # append new connection to temporary array
        end
      end
    end
    new_molecule.push(tempArray)   # add whole new atom connection list for atom to temporary array
  end
  new_molecule.each do |atom|   # again, sort connection numbers of each atom left to right from large to small
    elementSymbol=atom[0]
    atom.shift
    if(atom.length > 1)
      atom.sort!.reverse!
    end
    atom.insert(0,elementSymbol)
  end

  print "atomCount: ",atomCount,"\n\n"
  i=0
  exitLoop = false
  while(!exitLoop)
    print i,": ",new_molecule[i]," ",i+1,": ",new_molecule[i+1],"\n"
    if(new_molecule[i][0] == new_molecule[i+1][0])
      if(new_molecule[i].length > new_molecule[i+1].length)
        oldAtom=i
        newAtom=i+1
        print "Swap ",oldAtom," vs. ",newAtom,"\n\n"
        exitLoop = true
      end
    end
    i=i+1
    if(i >= atomCount)
      exitLoop = true
      print "\n"
      return new_molecule
    end
  end
  old_molecule=Marshal.load(Marshal.dump(new_molecule))
  old_molecule[oldAtom], old_molecule[newAtom] = old_molecule[newAtom], old_molecule[oldAtom]
  for i in 0..atomCount
    line=Array.new
    line=old_molecule[i]
    print line," "
    for j in 1..line.length-1
      print j,":",line[j]
      if(line[j] == oldAtom)
        print "c",newAtom
        line[j]=newAtom
      elsif(line[j] == newAtom)
        print "c",oldAtom
        line[j]=oldAtom
      end
      print " "
    end
    print "\n"
  end
  old_molecule.each do |atom|   # again, sort connection numbers of each atom left to right from large to small
    elementSymbol=atom[0]
    atom.shift
    if(atom.length > 1)
      atom.sort!.reverse!
    end
    atom.insert(0,elementSymbol)
  end
  print "\nSuggestion to re-order: ",old_molecule,"\n\n"
  return old_molecule
#
#
#
#    print "CorrespondenceTable [old,new]: \n\n",correspondenceTable,"\n\n"
#    print "New array WITH labels re-organized: \n\n",new_molecule,"\n\n"
#    return new_molecule
end


def canonicalization2(old_molecule)
  if(old_molecule == [])
    print "Structure is empty\n"
    return "Structure is empty"
  end
  periodicTable=Array.new
  periodicTable=PeriodicTable::Elements
  atomCount=old_molecule.length-1           # determine number of rows, make sure to properly account for array starting at index zero
  printf("Now sorting the array\n\n")
  old_molecule.each do |atom|   # first pass - sort connection numbers of each atom left to right from large to small
    elementSymbol=atom[0]
    atom.shift
    if(atom.length > 1)
      atom.sort!.reverse!
    end
    atom.insert(0,elementSymbol)
  end
  print "Old array: \n\n",old_molecule,"\n\n"
  new_molecule=Array.new
  correspondenceTable=Array.new   # contains pairs of array positions [old,new]
  j=0
  periodicTable.each do |element| # second pass - sort by element in increasing order, lowest atomic mass element to the left/bottom
    for i in 0..atomCount
      if(old_molecule[i][0] == element)
        new_molecule.push(old_molecule[i])
        correspondenceTable.push([i,j])
        j=j+1
      end
    end
  end
  old_molecule=Marshal.load(Marshal.dump(new_molecule))
  temp=Array.new
  temp=correspondenceTable.sort{|a,b| b[0] <=> a[0]}.reverse   # sort the correspondence table by lowest old atom position
  correspondenceTable=temp
  print "New array withOUT labels re-organized: \n\n",new_molecule,"\n\n"
  new_molecule=Array.new
  for i in 0..atomCount
    tempArray=Array.new
    tempArray.push(old_molecule[i][0])   # add element symbol to temporary array
    for j in 1..old_molecule[i].length-1
      for k in 0..correspondenceTable.length-1
        if(old_molecule[i][j] == correspondenceTable[k][0])
          tempArray.push(correspondenceTable[k][1])      # append new connection to temporary array
        end
      end
    end
    new_molecule.push(tempArray)   # add whole new atom connection list for atom to temporary array
  end
  new_molecule.each do |atom|   # again, sort connection numbers of each atom left to right from large to small
    elementSymbol=atom[0]
    atom.shift
    if(atom.length > 1)
      atom.sort!.reverse!
    end
    atom.insert(0,elementSymbol)
  end

  print "atomCount: ",atomCount,"\n\n"
  i=0
  exitLoop = false
  while(!exitLoop)
    print i,": ",new_molecule[i]," ",i+1,": ",new_molecule[i+1],"\n"
    if((new_molecule[i][0] == new_molecule[i+1][0]) && (new_molecule[i].length == new_molecule[i+1].length))
        do_swap=(new_molecule[i] <=> new_molecule[i+1])
#        do_swap=(new_molecule[i].to_s <=> new_molecule[i+1].to_s)
#        if(new_molecule[i].to_s > new_molecule[i+1].to_s)
      if(do_swap == 1)
        oldAtom=i
        newAtom=i+1
        print "Swap ",oldAtom," vs. ",newAtom,"\n\n"
        exitLoop = true
      end
    end
    i=i+1
    if(i >= atomCount)
      exitLoop = true
      print "\n"
      return new_molecule
    end
  end
  old_molecule=Marshal.load(Marshal.dump(new_molecule))
  old_molecule[oldAtom], old_molecule[newAtom] = old_molecule[newAtom], old_molecule[oldAtom]
  for i in 0..atomCount
    line=Array.new
    line=old_molecule[i]
    print line," "
    for j in 1..line.length-1
      print j,":",line[j]
      if(line[j] == oldAtom)
        print "c",newAtom
        line[j]=newAtom
      elsif(line[j] == newAtom)
        print "c",oldAtom
        line[j]=oldAtom
      end
      print " "
    end
    print "\n"
  end
  old_molecule.each do |atom|   # again, sort connection numbers of each atom left to right from large to small
    elementSymbol=atom[0]
    atom.shift
    if(atom.length > 1)
      atom.sort!.reverse!
    end
    atom.insert(0,elementSymbol)
  end
  print "\nSuggestion to re-order: ",old_molecule,"\n\n"
  return old_molecule
#
#
#
#    print "CorrespondenceTable [old,new]: \n\n",correspondenceTable,"\n\n"
#    print "New array WITH labels re-organized: \n\n",new_molecule,"\n\n"
#    return new_molecule
end


def serialization(molecule)
  if(molecule == [])
    print "\nStructure is empty\n"
    return "Structure is empty"
  elsif
    print molecule,"\n\n"
  end
  tempMolecule=Array.new
  tempMolecule=Marshal.load(Marshal.dump(molecule))
  graph=Array.new
  i=0
  tempMolecule.each do |atom|
    atom.shift
    atom.each do |connectedAtom|
      if(i < connectedAtom)
        graph.push([i,connectedAtom])
      else
        graph.push([connectedAtom,i])
      end
    end
    i=i+1
  end
  graph=graph.uniq.sort!
  print "Now printing new InChI\n\n"
  print "---------------------------------------\n"
  print "InChI=1S/sum_formula/"
  graph.each do |line|
    if(line[0] != line[1])
      print "(",line[0],"-",line[1],")"
    end
  end
  print "\n---------------------------------------\n"
  #
  # add output to inChIstring
  #
  inChIstring=""
  graph.each do |line|
    if(line[0] != line[1])
      inChIstring=inChIstring+"("+line[0].to_s+"-"+line[1].to_s+")"
    end
  end
  return inChIstring
end


def create_dot_file(molecule)
  dotfile=""
  print "\nPrinting Connection Table\n\n"
  if(molecule == [])
    print "Structure is empty\n"
    return "Structure is empty"
  elsif
    print molecule,"\n\n"
  end
  tempMolecule=Array.new
  tempMolecule=Marshal.load(Marshal.dump(molecule))
  graph=Array.new
  i=0
  tempMolecule.each do |atom|
    atom.shift
    atom.each do |connectedAtom|
      if(i < connectedAtom)
        graph.push([i,connectedAtom])
      else
        graph.push([connectedAtom,i])
      end
    end
    i=i+1
  end
  graph=graph.uniq.sort!
  elementColor=PeriodicTable::ElementColor
  #
  # print output to console
  #
  print "Now printing graph\n\n"
  print "---------------------------------------\n"
  print "graph test\n"
  print "{\n"
  print "  bgcolor=grey\n"
  i=0
  molecule.each do |atom|
    color = elementColor.fetch(:atom[0], "lightgrey")    # if color is unspecified fall back on lightgray
    print "  ",i," [label=\"",atom[0]," ",i,"\" color=",color,",style=filled,shape=circle,fontname=Calibri];\n"
    i=i+1
  end
  graph.each do |line|
    if(line[0] != line[1])
      print "  ",line[0]," -- ",line[1]," [color=black,style=bold];\n"
    end
  end
  print "}\n"
  print "---------------------------------------\n"
  #
  # add output to dotfile
  #
  dotfile=dotfile+"graph test\n"
  dotfile=dotfile+"{\n"
  dotfile=dotfile+"  bgcolor=grey\n"
  i=0
  molecule.each do |atom|
    color = elementColor.fetch(:atom[0], "lightgrey")
    dotfile=dotfile+"  "+i.to_s+" [label=\""+atom[0].to_s+" "+i.to_s+"\" color="+color+",style=filled,shape=circle,fontname=Calibri];\n"
    i=i+1
  end
  graph.each do |line|
    if(line[0] != line[1])
      dotfile=dotfile+"  "+line[0].to_s+" -- "+line[1].to_s+" [color=black,style=bold];\n"
    end
  end
  dotfile=dotfile+"}\n"
  return dotfile
end


def canonicalize_molecule(molecule)
  if(molecule.length > 1)
      for i in 0..(molecule.length-1)*20
          print "=== Start Pass #",i," ===\n\n"
          canonicalized_molecule = canonicalization1(molecule)
          canonicalized_molecule = canonicalization2(canonicalized_molecule)
          print "=== End Pass #",i," ===\n\n"
      end
  end
  return canonicalized_molecule
end


def create_ninchi_string(molecule)
  sumFormulaString=calculate_sum_formula(molecule)
  ninchi_string = "nInChI=1S/#{sumFormulaString}/c#{serialization(molecule)}"
  return ninchi_string
end
