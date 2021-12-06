#
# read input file and initialize data structures
#

def read_molfile(filename)
  molfile = File.read(filename) # reads entire file and closes it
  molfile.split("\n")
end

def create_node_features_matrix(molfile_lines, atom_count, periodic_table_elements)
  print "\nNow creating node_features_matrix ...\n"
  rows, columns, default_value = atom_count, 5, 0
  node_features_matrix = Array.new(rows) { Array.new(columns, default_value) }
  node_index = 0
  (7..atom_count + 6).each_with_index do |atom_index|
    line = molfile_lines[atom_index].split(' ')
    # node_features_matrix[node_index][0] is the atomic number
    # node_features_matrix[node_index][1] is the connectivity, which is added later in initialize_matrix()
    # node_features_matrix[node_index][2] is the "non-H" connectivity, which is added later in initialize_matrix()
    # node_features_matrix[node_index][3] is the "connectivity_index" which is added later in initialize_matrix() and then updated in sort_by_connectivity_index()
    # node_features_matrix[node_index][4] is the rest of the atom definition line, including the "pseudo-coordinates"
    # in general, further "node_features_matrix" elements starting from node_features_matrix[node_index][4] could be added here e.g. for isotopes, stereochemistry, ...
    node_features_matrix[node_index][0] = periodic_table_elements.index(line[3]) + 1
    n = line.length
    s = ''
    (4..n - 1).each do |i|
      s += " " + line[i]
    end
    node_features_matrix[node_index][4] = s
    node_index += 1
  end
  node_features_matrix
end

def create_edge_features_matrix(molfile_lines, edge_count, atom_count)
  print "\nNow creating edge_features_matrix ...\n"
  edge_features_matrix = Array.new(atom_count).map(&:to_a)
  (0..edge_count - 1).each do |edge_index|
    vertex1, vertex2 = parse_edge(molfile_lines[edge_index + 9 + atom_count])
    edge_features_matrix[vertex1].push(vertex2) # add to the first atom of a bond
    edge_features_matrix[vertex2].push(vertex1) # and to the second atom of the bond
  end
  edge_features_matrix
end

def parse_edge(molfile_line)
  vertex1 = molfile_line.split(' ')[4].to_i - 1 # need to substract 1 since molfile index starts at one, not at zero
  vertex2 = molfile_line.split(' ')[5].to_i - 1 # need to substract 1 since molfile index starts at one, not at zero
  vertex1, vertex2 = vertex2, vertex1 if vertex1 > vertex2 # make sure first atom always has lower (not: higher?) index
  [vertex1, vertex2]
end

def create_distance_matrix(adjacency_matrix)
  print "\nNow creating distance_matrix ...\n"
  n = adjacency_matrix.length
  rows, columns, default_value = n, n, 0
  distance_matrix = Array.new(rows) { Array.new(columns, default_value) }
  distance_matrix = adjacency_matrix
  if (n > 1)
    product_matrix = Array.new(rows) { Array.new(columns, default_value) }
    product_matrix = multiply_matrix(adjacency_matrix, adjacency_matrix)
    distance_matrix = update_distance_matrix(distance_matrix, product_matrix, 2)
  end
  if (n > 2)
    (3..n - 1).each do |i|
      product_matrix = multiply_matrix(adjacency_matrix, product_matrix)
      distance_matrix = update_distance_matrix(distance_matrix, product_matrix, i)
      print i, " "
    end
  end
  (0..n - 1).each do |i|
    distance_matrix[i][i] = 0 # diagonal elements should be zero, as entry is distance of atom to itself
  end
  distance_matrix
end

def multiply_matrix(matrix1, matrix2)
  n = matrix1.length
  rows, columns, default_value = n, n, 0
  product_matrix = Array.new(rows) { Array.new(columns, default_value) }
  (0..n - 1).each do |i|
    (0..n - 1).each do |j|
      (0..n - 1).each do |s|
        product_matrix[i][j] += matrix1[i][s] * matrix2[s][j]
      end
    end
  end
  product_matrix
end

def update_distance_matrix(old, new, d)
  n = old.length
  rows, columns, default_value = n, n, 0
  product_matrix = Array.new(rows) { Array.new(columns, default_value) }
  (0..n - 1).each do |i|
    (0..n - 1).each do |j|
      product_matrix[i][j] = old[i][j]
      if ((old[i][j] == 0) && (new[i][j] > 0))
        product_matrix[i][j] = d
      end
    end
  end
  product_matrix
end

def initialize_matrix(molfile_lines, periodic_table_elements)
  atom_count = molfile_lines[5].split(' ')[3].to_i # in molfile v3000, on 6th line, 1st number is number of atoms.
  edge_count = molfile_lines[5].split(' ')[4].to_i # in molfile v3000, on 6th line, 2nd number is number of bonds.
  rows, columns, default_value = atom_count, atom_count, 0
  adjacency_matrix = Array.new(rows) { Array.new(columns, default_value) }
  distance_matrix = Array.new(rows) { Array.new(columns, default_value) }
  edge_features_matrix = create_edge_features_matrix(molfile_lines, edge_count, atom_count)
  node_features_matrix = create_node_features_matrix(molfile_lines, atom_count, periodic_table_elements)
  molfile_header = Array.new(6) # the molfile v3000 header is six lines long
  (0..4).each do |line|
    molfile_header[line] = molfile_lines[line]
  end
  print "\nNow creating #{atom_count} x #{atom_count} adjacency matrix ...\n"
  (0..atom_count - 1).each do |row|
    line = edge_features_matrix[row]
    (0..atom_count - 1).each do |column|
      line.each do |entry|
        if (entry == column)
          adjacency_matrix[row][column] = 1 # here, in general, the edges (=bonds) could also be assigned additional properties by setting a value larger than 1, such as bond type/bond order (or a "bit field" or a another list/array)
          node_features_matrix[row][1] += 1 # count all neighbours here
        end
        if ((entry == column) && (node_features_matrix[column][0] > 1))
          node_features_matrix[row][2] += 1 # count only non-H neighbours here
        end
        node_features_matrix[row][3] =
          calculate_connectivity_index(adjacency_matrix, node_features_matrix, distance_matrix, row)
      end
    end
    print " ", row, " "
  end
  distance_matrix = create_distance_matrix(adjacency_matrix)
  offset = 6 + atom_count + 2 + edge_count + 2 # header is six lines, atom_block and bond_block enclosed in two tag lines each, which additionally have to be substracted
  line_count = molfile_lines.length - offset # header is six lines, atom_block and bond_block enclosed in two tag lines each, which additionally have to be substracted
  molfile_footer = Array.new(line_count)
  (0..line_count - 1).each do |line|
    molfile_footer[line] = molfile_lines[line + offset]
  end
  [adjacency_matrix, node_features_matrix, distance_matrix, molfile_header, molfile_footer]
end

#
# canonicalization and associated functions
#

def swap_matrix_elements(adjacency_matrix, node_features_matrix, distance_matrix, i, j)
  atom_count = node_features_matrix.length
  # swap rows of node_features_matrix
  atom_a = node_features_matrix[i]
  atom_b = node_features_matrix[j]
  node_features_matrix[i] = atom_b
  node_features_matrix[j] = atom_a
  # swap rows of adjacency_matrix
  for column in 0..atom_count-1
    atom_a = adjacency_matrix[i][column]
    atom_b = adjacency_matrix[j][column]
    adjacency_matrix[i][column] = atom_b
    adjacency_matrix[j][column] = atom_a
  end
  for row in 0..atom_count-1
    atom_a = adjacency_matrix[row][i]
    atom_b = adjacency_matrix[row][j]
    adjacency_matrix[row][i] = atom_b
    adjacency_matrix[row][j] = atom_a
  end
  # swap rows of distance_matrix
  for column in 0..atom_count-1
    atom_a = distance_matrix[i][column]
    atom_b = distance_matrix[j][column]
    distance_matrix[i][column] = atom_b
    distance_matrix[j][column] = atom_a
  end
  for row in 0..atom_count-1
    atom_a = distance_matrix[row][i]
    atom_b = distance_matrix[row][j]
    distance_matrix[row][i] = atom_b
    distance_matrix[row][j] = atom_a
  end
  [adjacency_matrix, node_features_matrix, distance_matrix]
end

def calculate_connectivity_index(adjacency_matrix, node_features_matrix, distance_matrix, row)
  atom_count = node_features_matrix.length
  connectivity_index = 0
  (0..atom_count - 1).each do |column|
    if (adjacency_matrix[row][column] > 0)
      connectivity_index += node_features_matrix[column][0] # the higher the atomic number of the neighbours the higher the assigned index number
    end
  end
  connectivity_index
end

def calculate_ac_index(node_features_matrix) # ac_index is a combination of atomic number and number of edges weighted by row number
  atom_count = node_features_matrix.length
  ac_index = 0
  for row in 0..(atom_count - 1)
    ac_index += node_features_matrix[row][0] * (row + 1) + node_features_matrix[row][2] * (row + 1)
  end
  ac_index
end

def sort_by_element_and_connectivity(adjacency_matrix, node_features_matrix, distance_matrix)
  print "\nNow sorting by atomic number and connectivity: \n"
  atom_count = node_features_matrix.length
  iteration = 1
  old_ac_index = 0
  new_ac_index = calculate_ac_index(node_features_matrix)
  while (old_ac_index < new_ac_index)
    for row in 0..(atom_count - 2)
      a = convert_to_bin(node_features_matrix[row][0], 8) + "-" + convert_to_bin(node_features_matrix[row][2], 5)
      b = convert_to_bin(node_features_matrix[row + 1][0],
                         8) + "-" + convert_to_bin(node_features_matrix[row + 1][2], 5)
      if (a > b)
        adjacency_matrix, node_features_matrix, distance_matrix = swap_matrix_elements(adjacency_matrix,
                                                                                       node_features_matrix, distance_matrix, row, row + 1)
      end
    end
    old_ac_index = new_ac_index
    new_ac_index = calculate_ac_index(node_features_matrix)
    print "\nIteration ", iteration, ": ", old_ac_index, " -> ", new_ac_index, "\n"
    iteration += 1
  end
  print "\n"
  print_adjacency_matrix(adjacency_matrix, node_features_matrix, distance_matrix)
  [adjacency_matrix, node_features_matrix, distance_matrix]
end

def sort_by_connectivity_index(adjacency_matrix, node_features_matrix, distance_matrix) # sort by connectivity index
  print "\nNow sorting by connectivity index: \n"
  iteration = 1
  converged = false
  atom_count = node_features_matrix.length
  previous_molecule_states = [Marshal.load(Marshal.dump(adjacency_matrix))]
  until converged == true
    print "\nIteration: #{iteration}\n"
    for row in 0..atom_count - 2
      for column in 0..atom_count - 1
        node_features_matrix[row][3] =
          calculate_connectivity_index(adjacency_matrix, node_features_matrix, distance_matrix, row)
        node_features_matrix[row + 1][3] =
          calculate_connectivity_index(adjacency_matrix, node_features_matrix, distance_matrix, row + 1)
      end
      if ((node_features_matrix[row][0] == node_features_matrix[row + 1][0]) && (node_features_matrix[row][1] == node_features_matrix[row + 1][1]) && (node_features_matrix[row][3] > node_features_matrix[row + 1][3]))
        adjacency_matrix, node_features_matrix, distance_matrix = swap_matrix_elements(adjacency_matrix,
                                                                                       node_features_matrix, distance_matrix, row, row + 1)
      end
    end
    if (previous_molecule_states.include?(adjacency_matrix))
      converged = true
    end
    previous_molecule_states.push(Marshal.load(Marshal.dump(adjacency_matrix)))
    iteration += 1
  end
  print "\n"
  print_adjacency_matrix(adjacency_matrix, node_features_matrix, distance_matrix)
  [adjacency_matrix, node_features_matrix, distance_matrix]
end

def sort_by_distance(adjacency_matrix, node_features_matrix, distance_matrix) # sort by distance to highest priority atom
  print "\nNow sorting by distance: \n\n"
  atom_count = node_features_matrix.length
  for i in 0..atom_count
    for j in 0..atom_count
      for row in 0..atom_count-2
        distance_A = distance_matrix[row][atom_count-1]
        distance_B = distance_matrix[row+1][atom_count-1]
        if((node_features_matrix[row][0] == node_features_matrix[row+1][0]) && (node_features_matrix[row][2] == node_features_matrix[row+1][2]) && (distance_A < distance_B))
          adjacency_matrix, node_features_matrix, distance_matrix = swap_matrix_elements(adjacency_matrix, node_features_matrix, distance_matrix, row, row+1)
          print "."
        end
      end
    end
  end
  print "\n"
  [adjacency_matrix, node_features_matrix, distance_matrix]
end

def sort_by_distance_index(adjacency_matrix, node_features_matrix, distance_matrix) # sort by distance to highest priority atom
  print "\nNow sorting by distance index: \n"
  atom_count = node_features_matrix.length
  for i in 0..atom_count
    for row in 0..atom_count-2
      for column in 0..atom_count-1
        distance_A = distance_matrix[row][atom_count-1]
        distance_B = distance_matrix[row+1][atom_count-1]
        distance_index_A = distance_matrix[row][column]*column
        distance_index_B = distance_matrix[row+1][column]*column
        if((node_features_matrix[row][0] == node_features_matrix[row+1][0]) && (node_features_matrix[row][2] == node_features_matrix[row+1][2]) && (distance_A == distance_B) && (distance_index_A > distance_index_B))
          adjacency_matrix, node_features_matrix, distance_matrix = swap_matrix_elements(adjacency_matrix, node_features_matrix, distance_matrix, row, row+1)
          print "."
        end
      end
    end
  end
  print "\n"
  [adjacency_matrix, node_features_matrix, distance_matrix]
end

def sort_terminal_hydrogens(adjacency_matrix, node_features_matrix, distance_matrix) # terminal hydrogen atom with highest index number to be attached to heavy atom with highest index number
  print "\nNow sorting terminal H atoms: \n"
  atom_count = node_features_matrix.length
  number_of_terminal_hydrogens = 0
  for row in 0..atom_count-1
    if((node_features_matrix[row][0] == 1) && (node_features_matrix[row][1] == 1))
      number_of_terminal_hydrogens += 1
    end
  end
  print "\nNumber of terminal H atoms: ",number_of_terminal_hydrogens,"\n"
  if(number_of_terminal_hydrogens < atom_count) # test required as otherwise program will crash for H2
    for i in 0..number_of_terminal_hydrogens-1
      print "\nIteration: ",i+1,"\n"
      for row in 0..number_of_terminal_hydrogens-2
        neighbour_A = 0
        neighbour_B = 0
        for column in 0..atom_count-1
          if(adjacency_matrix[row][column] != 0)
            neighbour_A = column
          end
          if(adjacency_matrix[row+1][column] != 0)
            neighbour_B = column
          end
        end
        if(neighbour_A > neighbour_B)
          adjacency_matrix, node_features_matrix, distance_matrix = swap_matrix_elements(adjacency_matrix, node_features_matrix, distance_matrix, row, row+1)
        end
      end
    end
  end
  [adjacency_matrix, node_features_matrix, distance_matrix]
end

def convert_to_bin(decimal_number, length_of_binary_number)
  binary_string = decimal_number.to_s(2)
  while (binary_string.length < length_of_binary_number)
    binary_string.insert(0, "0")
  end
  binary_string
end

def sort_adjacency_matrix(adjacency_matrix, node_features_matrix, distance_matrix)
  atom_count = node_features_matrix.length
  print "\nNow sorting adjacency matrix: \n\n"
  print_adjacency_matrix(adjacency_matrix, node_features_matrix, distance_matrix)
  adjacency_matrix, node_features_matrix, distance_matrix = sort_by_element_and_connectivity(adjacency_matrix, node_features_matrix, distance_matrix)
  adjacency_matrix, node_features_matrix, distance_matrix = sort_by_connectivity_index(adjacency_matrix, node_features_matrix, distance_matrix)
  adjacency_matrix, node_features_matrix, distance_matrix = sort_by_distance(adjacency_matrix, node_features_matrix, distance_matrix)
  adjacency_matrix, node_features_matrix, distance_matrix = sort_by_distance_index(adjacency_matrix, node_features_matrix, distance_matrix)
  adjacency_matrix, node_features_matrix, distance_matrix = sort_terminal_hydrogens(adjacency_matrix, node_features_matrix, distance_matrix)  
  [adjacency_matrix, node_features_matrix, distance_matrix]
end

#
# serialization and various other output, including *.dot and *.mol files
#

def print_adjacency_matrix(adjacency_matrix, node_features_matrix, distance_matrix)
  atom_count = adjacency_matrix.length
  distance_index = 0
  (0..atom_count - 1).each do |row|
    print adjacency_matrix[row], " ", distance_matrix[row], " ", node_features_matrix[row][0..3]
    neighbours = ''
    (0..atom_count - 1).each do |column|
      if (adjacency_matrix[row][column] == 1)
        neighbours = neighbours + column.to_s + ","
      end
      distance_index += distance_matrix[row][column] * column
    end
    neighbours.chop!
    print " {#{neighbours}} {#{distance_index}} " # {}"\n"
    # distance_index = 0
    print convert_to_bin(node_features_matrix[row][0], 8), " - ", convert_to_bin(node_features_matrix[row][2], 5), " - ",
          convert_to_bin(distance_index, distance_matrix.length), " \n"
    distance_index = 0
  end
end

def write_ninchi_string(adjacency_matrix, node_features_matrix, periodic_table_elements)
  sum_formula = write_sum_formula_string(node_features_matrix, periodic_table_elements)
  serialized_molecule = serialize_molecule(adjacency_matrix, node_features_matrix)
  "nInChI=1S/#{sum_formula}/c#{serialized_molecule}"
end

def write_dot_file(adjacency_matrix, node_features_matrix, filename, periodic_table_elements, periodic_table_colors)
  filename = File.basename(filename, '.mol')
  dotfile = "graph #{filename}\n{\n  bgcolor=grey\n"
  n = adjacency_matrix.length
  (0..n - 1).each do |i|
    symbol = periodic_table_elements[node_features_matrix[i][0] - 1]
    color = periodic_table_colors.fetch(symbol, 'lightgrey')
    dotfile += "  #{i} [label=\"#{symbol} #{i}\" color=#{color},style=filled,shape=circle,fontname=Calibri];\n"
  end
  (0..n - 1).each do |row|
    (0..n - 1).each do |column|
      if (adjacency_matrix[row][column] == 1)
        if (row < column)
          dotfile += "  #{row.to_s} -- #{column.to_s} [color=black,style=bold];\n"
        end
      end
    end
  end
  dotfile += "}\n"
end

def write_molfile(adjacency_matrix, node_features_matrix, molfile_header, molfile_footer, periodic_table_elements)
  molfile = ''
  # print header block
  molfile_header[1] = "nInChI v2.5"
  (0..4).each do |line|
    molfile += molfile_header[line] + "\n"
  end
  atom_count = node_features_matrix.length
  bond_count = 0
  (0..atom_count - 1).each do |row|
    bond_count += node_features_matrix[row][1].to_i
  end
  bond_count = bond_count / 2 # divide by two since each bond is present twice in the node_features_matrix, for both A and B of an A-B bond
  molfile += "M  V30 COUNTS " + atom_count.to_s + " " + bond_count.to_s + " 0 0 1\n"
  # print atom block
  molfile += "M  V30 BEGIN ATOM\n"
  (0..atom_count - 1).each do |i|
    symbol = periodic_table_elements[node_features_matrix[i][0] - 1]
    molfile += "M  V30 " + (i + 1).to_s + " " + symbol + node_features_matrix[i][4] + "\n"
  end
  molfile += "M  V30 END ATOM\n"
  molfile += "M  V30 BEGIN BOND\n"
  # print bond block
  bond_count = 0
  bond_order = 1
  (0..atom_count - 1).each do |row|
    (0..atom_count - 1).each do |column|
      if ((row < column) && (adjacency_matrix[row][column] == 1))
        molfile += "M  V30 " + (bond_count + 1).to_s + " " + bond_order.to_s + " " + (row + 1).to_s + " " + (column + 1).to_s + " 0 0 0 0\n"
        bond_count += 1
      end
    end
  end
  molfile += "M  V30 END BOND"
  # print rest of molfile
  molfile_footer.each do |line|
    molfile += "\n" + line
  end
  molfile += "\n"
  molfile += "\nWARNING: atom index numbers beyond bond block are currently not updated!!!"
end

def standard_inchi(adjacency_matrix, node_features_matrix)
  inchi_string = '/c'
  inchi_string_H = ''
  n = adjacency_matrix.length
  (0..n - 1).each do |row|
    (0..n - 1).each do |column|
      if ((adjacency_matrix[row][column] == 1) && (row < column))
        if (node_features_matrix[row][0] != 1)
          inchi_string += '(' + row.to_s + '-' + column.to_s + ')'
        end
        if (node_features_matrix[row][0] == 1)
          inchi_string_H += '(' + row.to_s + '-' + column.to_s + ')'
        end
      end
    end
  end
  if (inchi_string_H != '')
    inchi_string = inchi_string + '/h' + inchi_string_H
  end
  inchi_string
end

def serialize_molecule(adjacency_matrix, node_features_matrix)
  #
  # inchi_string[0] -> "2-tuple format" (0-2)(1-2)
  # inchi_string[1] -> "n-tuple format" (2:0,1)
  # inchi_string[2] -> binary format with each row of adjacency matrix elements separated by minus sign
  # inchi_string[3] -> binary number generated from inchi_string[2] by removal of all "-" and addition of leading "1"
  # inchi_string[4] -> decimal format of binary inchi_string[3]
  # inchi_string[5] -> hexadecimal format of binary inchi_string[3]
  # inchi_string[6] -> base32 encoding of binary inchi_string[3]
  # inchi_string[7] -> attempt to "InChI-sytle" separate output of heavy vs. hydrogen atoms
  # inchi_string[8] -> sum of neighbour atom index numbers to detect degeneracy problem
  #
  inchi_string = Array.new(9, '')
  n = adjacency_matrix.length
  (0..n - 1).each do |row|
    inchi_string[1] += '(' + row.to_s + ':'
    inchi_string[8] += '(' + node_features_matrix[row][0].to_s + "_" + row.to_s + ':'
    index_sum = 0
    (0..n - 1).each do |column|
      inchi_string[2] += adjacency_matrix[row][column].to_s
      if (adjacency_matrix[row][column] == 1)
        inchi_string[1] += column.to_s + ','
        index_sum = index_sum + column
        if (row < column)
          inchi_string[0] += '(' + row.to_s + '-' + column.to_s + ')'
        end
      end
    end
    inchi_string[8] += index_sum.to_s
    inchi_string[2] += '-'
    inchi_string[1].chop!
    inchi_string[1] += ')'
    inchi_string[8] += ')'
  end
  inchi_string[2].chop!
  inchi_string[3] += '1' + inchi_string[2] # add a leading "1" to avoid numbers with leading zeroes to become equivalent
  inchi_string[3].delete! '-'
  inchi_string[4] = inchi_string[3].to_i(2)
  inchi_string[5] = 'hex:' + inchi_string[4].to_s(base = 16)
  inchi_string[6] = 'base32:' + inchi_string[4].to_s(base = 32)
  inchi_string[7] = standard_inchi(adjacency_matrix, node_features_matrix)
  inchi_string[0] # only return tuples for now (makes reading output more informative)
end

def write_sum_formula_string(node_features_matrix, periodic_table_elements)
  # Write sum formula in the order C > H > all other elements in alphabetic order.
  element_counts = compute_element_counts(node_features_matrix, periodic_table_elements)
  element_counts.transform_values! { |v| v > 1 ? v : '' } # remove 1s since counts of 1 are implicit in sum formula
  sum_formula_string = ''
  sum_formula_string += "C#{element_counts['C']}" if element_counts.key?('C')
  sum_formula_string += "H#{element_counts['H']}" if element_counts.key?('H')
  element_counts.sort.to_h.each do |element, count|
    sum_formula_string += "#{element}#{count}" unless %w[C H].include?(element)
  end
  sum_formula_string
end

def compute_element_counts(node_features_matrix, periodic_table_elements)
  # Compute hash table mapping element symbols to stoichiometric counts.
  atom_list = []
  node_features_matrix.each { |atom| atom_list.push(atom[0]) }
  unique_elements = atom_list.uniq
  initial_counts = Array.new(unique_elements.length, 0)
  element_counts = unique_elements.zip(initial_counts).to_h
  atom_list.each { |atom| element_counts[atom] += 1 }
  element_counts.transform_keys! { |k|
    periodic_table_elements[k - 1]
  } # change atomic mass to element symbol; k - 1; convert from mass back to index
end
