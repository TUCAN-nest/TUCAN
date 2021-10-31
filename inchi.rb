# (c) CC BY-SA | Jan C. Brammer, RWTH Aachen and Ulrich Schatzschneider, Universität Würzburg | NFDI4Chem | v2.1 | 18.09.2021

require './periodic_table'

module Inchi
  def read_molfile(filename)
    molfile = File.read(filename) # reads entire file and closes it
    molfile.split("\n")
  end

  def create_node_features_matrix(molfile_lines, atom_count, periodic_table_elements)
    rows, columns, default_value = atom_count, 3, 0
    node_features_matrix = Array.new(rows) { Array.new(columns, default_value) }
    node_index = 0
    (7..atom_count + 6).each_with_index do |atom_index|
      atom = molfile_lines[atom_index].split(' ')[3]
      # node_features_matrix.push(periodic_table_elements.index(atom) + 1)
      node_features_matrix[node_index][0] = periodic_table_elements.index(atom) + 1
      node_index += 1
    end
    node_features_matrix
  end

  def create_edge_features_matrix(molfile_lines, edge_count, atom_count)
    edge_features_matrix = Array.new(atom_count).map(&:to_a)
    (0..edge_count - 1).each do |edge_index|
      vertex1, vertex2 = parse_edge(molfile_lines[edge_index + 9 + atom_count])
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

  def create_adjacency_matrix(molfile_lines, periodic_table_elements)
    atom_count = molfile_lines[5].split(' ')[3].to_i # in molfile v3000, on 6th line, 1st number is number of atoms.
    edge_count = molfile_lines[5].split(' ')[4].to_i # in molfile v3000, on 6th line, 2nd number is number of bonds.
    node_features_matrix = create_node_features_matrix(molfile_lines, atom_count, periodic_table_elements)
    edge_features_matrix = create_edge_features_matrix(molfile_lines, edge_count, atom_count)
    rows, columns, default_value = atom_count, atom_count, 0
    adjacency_matrix = Array.new(rows) { Array.new(columns, default_value) }
    # in general, further "atom property lists" could be added here e.g. for isotopes, stereochemistry, ...
    print "\nAdjacency matrix is #{atom_count} x #{atom_count}\n"
    (0..atom_count - 1).each do |row|
      line = edge_features_matrix[row]
      (0..atom_count - 1).each do |column|
        line.each do |entry|
          if (entry == column)
            adjacency_matrix[row][column] = 1 # here, in general, the edges (=bonds) could also be assigned additional properties by setting a value larger than 1, such as bond type/bond order (or a "bit field" or a another list/array)
            node_features_matrix[row][1] += 1
          end
        end
      end
    end
    [adjacency_matrix, node_features_matrix] # better also directly calculate and return "neighbour list" here, as needed in the 3rd sorting step
  end

  def swap_adjacency_matrix_elements(adjacency_matrix, node_features_matrix, i, j)
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
    [adjacency_matrix, node_features_matrix]
  end

def sort_by_element(adjacency_matrix, node_features_matrix) # sort by atomic mass
    atom_count = node_features_matrix.length
    for i in (atom_count).downto(2)
      for j in 0..(atom_count-2)
        if(node_features_matrix[j][0] > node_features_matrix[j+1][0])
          adjacency_matrix, node_features_matrix = swap_adjacency_matrix_elements(adjacency_matrix, node_features_matrix, j, j+1)
        end
      end
    end
    [adjacency_matrix, node_features_matrix]
end

def sort_by_connectivity(adjacency_matrix, node_features_matrix) # sort by connectivity
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
    [adjacency_matrix, node_features_matrix]
end

def sort_by_connectivity_index(adjacency_matrix, node_features_matrix) # sort by connectivity index
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
          if(adjacency_matrix[row][column] == 1)
            connectivity_index_A = connectivity_index_A+adjacency_matrix[row][column]*(column+1) # column index number * matrix element (0 or 1)
          end
          if(adjacency_matrix[row+1][column] == 1)
            connectivity_index_B = connectivity_index_B+adjacency_matrix[row+1][column]*(column+1) # column index number * matrix element (0 or 1)
          end
        end
        if((node_features_matrix[row][0] == node_features_matrix[row+1][0]) && (node_features_matrix[row][1] == node_features_matrix[row+1][1]) && (connectivity_index_B < connectivity_index_A))
          adjacency_matrix, node_features_matrix = swap_adjacency_matrix_elements(adjacency_matrix, node_features_matrix, row, row+1)
        end
      end

      if(previous_molecule_states.include?(adjacency_matrix))
        converged = true
        print "\nOptimization has converged\n"
      end

      previous_molecule_states.push(Marshal.load(Marshal.dump(adjacency_matrix)))

      iteration += 1
      print_adjacency_matrix(adjacency_matrix, node_features_matrix)
    end
    [adjacency_matrix, node_features_matrix]
  end

  def sort_adjacency_matrix(adjacency_matrix, node_features_matrix)
    atom_count = node_features_matrix.length
    print "\nNow sorting adjacency matrix\n"
    print_adjacency_matrix(adjacency_matrix, node_features_matrix)
    print "\nNumber of atoms: #{atom_count}\n"
    adjacency_matrix, node_features_matrix = sort_by_element(adjacency_matrix, node_features_matrix)
    adjacency_matrix, node_features_matrix = sort_by_connectivity(adjacency_matrix, node_features_matrix)
    adjacency_matrix, node_features_matrix = sort_by_connectivity_index(adjacency_matrix, node_features_matrix)
    [adjacency_matrix, node_features_matrix]
  end

  def print_adjacency_matrix(adjacency_matrix, node_features_matrix)
    n = adjacency_matrix.length
    (0..n - 1).each do |row|
      print " ",adjacency_matrix[row]," ",node_features_matrix[row]
      connectivity_index = 0 # need to set back to zero for each new row
      neighbours = ''
      print "\n["
      (0..n - 1).each do |column|
        connectivity_index = connectivity_index + adjacency_matrix[row][column] * (column + 1) # column index number * matrix element (0 or 1)
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
    #
    inchi_string = Array.new(8, '')
    n = adjacency_matrix.length
    (0..n - 1).each do |row|
      inchi_string[1] += '(' + row.to_s + ':'
      (0..n - 1).each do |column|
        inchi_string[2] += adjacency_matrix[row][column].to_s
        if (adjacency_matrix[row][column] == 1)
          inchi_string[1] += column.to_s + ','
          if (row < column)
            inchi_string[0] += '(' + row.to_s + '-' + column.to_s + ')'
          end
        end
      end
      inchi_string[2] += '-'
      inchi_string[1].chop!
      inchi_string[1] += ')'
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
