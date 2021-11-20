# (c) CC BY-SA | Jan C. Brammer, RWTH Aachen and Ulrich Schatzschneider, Universität Würzburg | NFDI4Chem | v2.5 | 20.11.2021

require './periodic_table'

module Inchi
  def read_molfile(filename)
    molfile_lines = File.read(filename).split("\n") # reads entire file and closes it
    atom_count = molfile_lines[5].split(' ')[3].to_i # in molfile v3000, on 6th line, 1st number is number of atoms.
    edge_count = molfile_lines[5].split(' ')[4].to_i # in molfile v3000, on 6th line, 2nd number is number of bonds.
    atom_block = (7..atom_count + 6).map { |i| molfile_lines[i] }
    edge_block = (0..edge_count - 1).map { |i| molfile_lines[i + 9 + atom_count] }
    molfile_header = Array.new(6) # the molfile v3000 header is six lines long
    (0..4).each do |line|
      molfile_header[line] = molfile_lines[line]
    end
    [atom_block, edge_block, molfile_lines, molfile_header]
  end

  def create_node_features_matrix(atom_block, periodic_table_elements)
    rows, columns, default_value = atom_block.size, 3, 0
    node_features_matrix = Array.new(rows) { Array.new(columns, default_value) }
    node_index = 0
    atom_block.each do |atom_line|
      symbol = atom_line.split(' ')[3]
      # node_features_matrix.push(periodic_table_elements.index(atom) + 1)
      node_features_matrix[node_index][0] = periodic_table_elements.index(symbol) + 1 # atomic number (German: Ordnungszahl)
      # node_features_matrix[node_index][1] is the connectivity, which is added later in create_adjacency_matrix()
      # node_features_matrix[node_index][2] is the "connectivity_index" which is also added later in create_adjacency_matrix() and then updated in sort_by_connectivity_index()
      # in general, further "node_features_matrix" elements starting from node_features_matrix[node_index][3] could be added here e.g. for isotopes, stereochemistry, ...
      node_index += 1
    end
    node_features_matrix
  end

  def create_edge_features_matrix(edge_block, atom_count)
    edge_features_matrix = Array.new(atom_count).map(&:to_a)
    edge_block.each do |edge_line|
      vertex1, vertex2 = parse_edge(edge_line)
      edge_features_matrix[vertex1].push(vertex2)    # add to the first atom of a bond
      edge_features_matrix[vertex2].push(vertex1)    # and to the second atom of the bond
    end
    edge_features_matrix
  end

  def parse_edge(molfile_line)
    vertex1 = molfile_line.split(' ')[4].to_i - 1 # need to substract 1 since molfile index starts at one, not at zero
    vertex2 = molfile_line.split(' ')[5].to_i - 1 # need to substract 1 since molfile index starts at one, not at zero
    vertex1, vertex2 = vertex2, vertex1 if vertex1 > vertex2 # make sure first atom always has lower (not: higher?) index
    [vertex1, vertex2]
  end

  def calculate_connectivity_index(adjacency_matrix, node_features_matrix, atom_count, row)
    connectivity_index = 0
    for column in 0..atom_count-1
      connectivity_index = connectivity_index+adjacency_matrix[row][column]*(column+1)*node_features_matrix[column][0]
    end
    connectivity_index
  end
  
  def initialize_matrix(atom_block, edge_block, header_block, periodic_table_elements)
    atom_count = atom_block.size
    rows, columns, default_value = atom_count, atom_count, 0
    adjacency_matrix = Array.new(rows) { Array.new(columns, default_value) }
    distance_matrix = Array.new(rows) { Array.new(columns, default_value) }
    edge_features_matrix = create_edge_features_matrix(edge_block, atom_count)
    node_features_matrix = create_node_features_matrix(atom_block, periodic_table_elements)
    molfile_header = Array.new(6) # the molfile v3000 header is six lines long
    (0..4).each do |line|
      molfile_header[line] = header_block[line]
    end
    print "\nAdjacency matrix is #{atom_count} x #{atom_count}\n"
    (0..atom_count - 1).each do |row|
      line = edge_features_matrix[row]
      (0..atom_count - 1).each do |column|
        line.each do |entry|
          if (entry == column)
            adjacency_matrix[row][column] = 1 # here, in general, the edges (=bonds) could also be assigned additional properties by setting a value larger than 1, such as bond type/bond order (or a "bit field" or a another list/array)
            node_features_matrix[row][1] += 1
          end
          node_features_matrix[row][2] = calculate_connectivity_index(adjacency_matrix, node_features_matrix, atom_count, row)
        end
      end
    end
    [adjacency_matrix, node_features_matrix, distance_matrix, molfile_header]
  end

  def swap_adjacency_matrix_elements(adjacency_matrix, node_features_matrix, distance_matrix, i, j)
    atom_a = node_features_matrix[i]
    atom_b = node_features_matrix[j]
    node_features_matrix[i] = atom_b
    node_features_matrix[j] = atom_a
    atom_count = node_features_matrix.length - 1
    for column in 0..atom_count
      atom_a = adjacency_matrix[i][column]
      atom_b = adjacency_matrix[j][column]
      adjacency_matrix[i][column] = atom_b
      adjacency_matrix[j][column] = atom_a
    end
    for row in 0..atom_count
      atom_a = adjacency_matrix[row][i]
      atom_b = adjacency_matrix[row][j]
      adjacency_matrix[row][i] = atom_b
      adjacency_matrix[row][j] = atom_a
    end
    [adjacency_matrix, node_features_matrix, distance_matrix]
  end

def sort_by_element(adjacency_matrix, node_features_matrix, distance_matrix) # sort by atomic mass
    atom_count = node_features_matrix.length
    for i in (atom_count).downto(2)
      for j in 0..(atom_count-2)
        if(node_features_matrix[j][0] > node_features_matrix[j+1][0])
          adjacency_matrix, node_features_matrix = swap_adjacency_matrix_elements(adjacency_matrix, node_features_matrix, j, j+1)
        end
      end
    end
    [adjacency_matrix, node_features_matrix, distance_matrix]
end

def sort_by_connectivity(adjacency_matrix, node_features_matrix, distance_matrix) # sort by connectivity
    atom_count = node_features_matrix.length
    for i in (atom_count).downto(2)
      for j in 0..(atom_count-2)
        for row in 0..atom_count-2
          if((node_features_matrix[row][0] == node_features_matrix[row+1][0]) && (node_features_matrix[row][1] > node_features_matrix[row+1][1]))
            adjacency_matrix, node_features_matrix = swap_adjacency_matrix_elements(adjacency_matrix, node_features_matrix, row, row+1)
          end
        end
      end
    end
    [adjacency_matrix, node_features_matrix, distance_matrix]
end

def sort_by_connectivity_index(adjacency_matrix, node_features_matrix, distance_matrix) # sort by connectivity index
    iteration = 1
    converged = false
    atom_count = node_features_matrix.length
    previous_molecule_states = [Marshal.load(Marshal.dump(adjacency_matrix))]
    until converged == true
      print "\nCycle ##{iteration}\n"
      for row in 0..atom_count-2
        connectivity_index_A = 0
        connectivity_index_B = 0
        for column in 0..atom_count-1
          connectivity_index_A = connectivity_index_A+adjacency_matrix[row][column]*(column+1)*node_features_matrix[column][0]
          connectivity_index_B = connectivity_index_B+adjacency_matrix[row+1][column]*(column+1)*node_features_matrix[column][0]
        end
        if((node_features_matrix[row][0] == node_features_matrix[row+1][0]) && (node_features_matrix[row][1] == node_features_matrix[row+1][1]) && (connectivity_index_B < connectivity_index_A))
          adjacency_matrix, node_features_matrix = swap_adjacency_matrix_elements(adjacency_matrix, node_features_matrix, row, row+1)
        end
        node_features_matrix[row][2] = calculate_connectivity_index(adjacency_matrix, node_features_matrix, atom_count, row)
        node_features_matrix[row+1][2] = calculate_connectivity_index(adjacency_matrix, node_features_matrix, atom_count, row+1)
      end

      if(previous_molecule_states.include?(adjacency_matrix))
        converged = true
        print "\nOptimization has converged\n"
      end

      previous_molecule_states.push(Marshal.load(Marshal.dump(adjacency_matrix)))

      iteration += 1
      print_adjacency_matrix(adjacency_matrix, node_features_matrix)
    end
    [adjacency_matrix, node_features_matrix, distance_matrix]
  end

  def sort_adjacency_matrix(adjacency_matrix, node_features_matrix, distance_matrix)
    atom_count = node_features_matrix.length
    print "\nNow sorting adjacency matrix\n"
    print_adjacency_matrix(adjacency_matrix, node_features_matrix)
    print "\nNumber of atoms: #{atom_count}\n"
    adjacency_matrix, node_features_matrix = sort_by_element(adjacency_matrix, node_features_matrix)
    adjacency_matrix, node_features_matrix = sort_by_connectivity(adjacency_matrix, node_features_matrix)
    adjacency_matrix, node_features_matrix = sort_by_connectivity_index(adjacency_matrix, node_features_matrix)
    [adjacency_matrix, node_features_matrix, distance_matrix]
  end

  def print_adjacency_matrix(adjacency_matrix, node_features_matrix, distance_matrix)
    n = adjacency_matrix.length
    (0..n - 1).each do |row|
      print " ",adjacency_matrix[row]," ",node_features_matrix[row]
      connectivity_index = 0 # need to set back to zero for each new row
      neighbours = ''
      (0..n - 1).each do |column|
        connectivity_index = connectivity_index + adjacency_matrix[row][column] * (column + 1) * node_features_matrix[column][0]
        if (adjacency_matrix[row][column] == 1)
          neighbours = neighbours + column.to_s + ","
        end
      end
      neighbours.chop!
      print " {#{connectivity_index}} {#{neighbours}}\n"
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
  
  def write_molfile(adjacency_matrix, node_features_matrix, molfile_header, periodic_table_elements)
    molfile = ''
    molfile_header[1] = "nInChI v2.5"
    (0..4).each do |line|
      molfile += molfile_header[line]+"\n"
    end
    atom_count = node_features_matrix.length
    bond_count = 0
    (0..atom_count - 1).each do |row|
      bond_count += node_features_matrix[row][1].to_i
    end
    bond_count = bond_count/2 # divide by two since each bond is present twice in the node_features_matrix, for both A and B of an A-B bond
    molfile += "M  V30 COUNTS "+atom_count.to_s+" "+bond_count.to_s+" 0 0 1\n"
    molfile += "M  V30 BEGIN ATOM\n"
    (0..atom_count-1).each do |i|
      symbol = periodic_table_elements[node_features_matrix[i][0]-1]
      molfile += "M  V30 "+(i+1).to_s+" "+symbol+" 0.0 0.0 0.0 0\n"
    end
    molfile += "M  V30 END ATOM\n"
    molfile += "M  V30 BEGIN BOND\n"
    # print bond block
    bond_count = 0
    bond_order = 1
    (0..atom_count - 1).each do |row|
      (0..atom_count - 1).each do |column|
        if ((row < column) && (adjacency_matrix[row][column] == 1))
          molfile += "M  V30 "+(bond_count+1).to_s+" "+bond_order.to_s+" "+(row+1).to_s+" "+(column+1).to_s+" 0 0 0 0\n"
          bond_count += 1
        end
      end
    end
    molfile += "M  V30 END BOND\n"
    molfile += "M  V30 END CTAB\n"
    molfile += "M  END\n"
  end
  
  def standard_inchi(adjacency_matrix, node_features_matrix)
    inchi_string = '/c'
    inchi_string_H = ''
    n = adjacency_matrix.length
    (0..n-1).each do |row|
     (0..n-1).each do |column|
        if((adjacency_matrix[row][column] == 1) && (row < column))
          if(node_features_matrix[row][0] != 1)
            inchi_string += '('+row.to_s+'-'+column.to_s+')'
          end
          if(node_features_matrix[row][0] == 1)
            inchi_string_H += '('+row.to_s+'-'+column.to_s+')'
          end
        end
      end
    end
    if(inchi_string_H != '')
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
      inchi_string[8] += '('+node_features_matrix[row][0].to_s+"_"+row.to_s+':'
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
    inchi_string[7]=standard_inchi(adjacency_matrix, node_features_matrix)
    inchi_string
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
    element_counts.transform_keys! { |k| periodic_table_elements[k - 1] } # change atomic mass to element symbol; k - 1; convert from mass back to index
  end
end
