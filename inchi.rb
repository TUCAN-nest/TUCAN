# (c) CC BY-SA | Jan C. Brammer, RWTH Aachen and Ulrich Schatzschneider, Universität Würzburg | NFDI4Chem | v2.1 | 18.09.2021

require './periodic_table'

module Inchi
  def read_molfile(filename)
    molfile = File.read(filename) # reads entire file and closes it
    molfile.split("\n")
  end

  def create_atom_list(molfile_lines, atom_count, periodic_table_elements)
    elements = []
    (4..atom_count + 3).each_with_index do |atom_index|
      atom = molfile_lines[atom_index].split(' ')[3]
      elements.push(periodic_table_elements.index(atom) + 1)
    end
    elements
  end

  def create_edge_list(molfile_lines, edge_count, atom_count)
    edges = Array.new(atom_count).map(&:to_a)
    (0..edge_count - 1).each do |edge_index|
      vertex1, vertex2 = parse_edge(molfile_lines[edge_index + 4 + atom_count])
      edges[vertex1].push(vertex2)    # add to the first atom of a bond
      edges[vertex2].push(vertex1)    # and to the second atom of the bond
    end
    edges
  end

  def parse_edge(molfile_line)
    vertex1, vertex2, * = molfile_line.split(' ').map { |i| i.to_i - 1 }
    vertex1, vertex2 = vertex2, vertex1 if vertex1 > vertex2 # make sure first atom always has lower (not: higher?) index
    [vertex1, vertex2]
  end

  def create_adjacency_matrix(molfile_lines, periodic_table_elements)
    atom_count, edge_count = molfile_lines[3].scan(/\d+/).map(&:to_i) # on 4th  line, 1st number is number of atoms, 2nd number is number of bonds.
    atom_list = create_atom_list(molfile_lines, atom_count, periodic_table_elements)
    edge_list = create_edge_list(molfile_lines, edge_count, atom_count)
    rows, columns, default_value = atom_count, atom_count, 0
    adjacency_matrix = Array.new(rows) { Array.new(columns, default_value) }
    # in general, further "atom property lists" could be added here e.g. for isotopes, stereochemistry, ...
    print "\nAdjacency matrix is #{atom_count} x #{atom_count}\n"
    (0..atom_count - 1).each do |row|
      line = edge_list[row]
      (0..atom_count - 1).each do |column|
        line.each do |entry|
          if (entry == column)
            adjacency_matrix[row][column] = 1 # here, in general, the edges (=bonds) could also be assigned additional properties by setting a value larger than 1, such as bond type/bond order (or a "bit field" or a another list/array)
          end
        end
      end
    end
    [adjacency_matrix, atom_list] # better also directly calculate and return "neighbour list" here, as needed in the 3rd sorting step
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
        if(node_features_matrix[j] > node_features_matrix[j+1])
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
          edge_count_A = 0
          edge_count_B = 0
          for column in 0..atom_count-1
            if(adjacency_matrix[row][column] == 1)
              edge_count_A += 1
            end
            if(adjacency_matrix[row+1][column] == 1)
              edge_count_B += 1
            end
          end
          if((node_features_matrix[row] == node_features_matrix[row+1]) && (edge_count_A > edge_count_B))
            adjacency_matrix, node_features_matrix = swap_adjacency_matrix_elements(adjacency_matrix, node_features_matrix, row, row+1)
          end
        end
      end
    end
    [adjacency_matrix, node_features_matrix]
end
  
  def sort_adjacency_matrix(adjacency_matrix, node_features_matrix)
    print "\nNow sorting adjacency matrix\n"
    print_adjacency_matrix(adjacency_matrix, node_features_matrix)
    atom_count = node_features_matrix.length
    print "\nNumber of atoms: #{atom_count}\n"
    iteration = 1
    converged = false
    adjacency_matrix, node_features_matrix = sort_by_element(adjacency_matrix, node_features_matrix)
    adjacency_matrix, node_features_matrix = sort_by_connectivity(adjacency_matrix, node_features_matrix)
    previous_molecule_states = [Marshal.load(Marshal.dump(adjacency_matrix))]
    until converged == true
      print "\nIteration ##{iteration}\n"
  
      #
      # sort by connectivity_index
      #

      for row in 0..atom_count - 2
        connectivity_index_A = 0
        connectivity_index_B = 0
        edge_count_A = 0
        edge_count_B = 0
        for column in 0..atom_count - 1
          if (adjacency_matrix[row][column] == 1)
            edge_count_A += 1
            connectivity_index_A = connectivity_index_A + adjacency_matrix[row][column] * (column + 1) # column index number * matrix element (0 or 1)
          end
          if (adjacency_matrix[row + 1][column] == 1)
            edge_count_B += 1
            connectivity_index_B = connectivity_index_B + adjacency_matrix[row + 1][column] * (column + 1) # column index number * matrix element (0 or 1)
          end
        end
        if ((node_features_matrix[row] == node_features_matrix[row + 1]) && (edge_count_A == edge_count_B) && (connectivity_index_B < connectivity_index_A))
          adjacency_matrix, node_features_matrix = swap_adjacency_matrix_elements(adjacency_matrix, node_features_matrix, row, row + 1)
        end
      end

      if (previous_molecule_states.include?(adjacency_matrix))
        converged = true
        print "\nSorting has converged.\n"
      end

      previous_molecule_states.push(Marshal.load(Marshal.dump(adjacency_matrix)))

      iteration += 1
      print_adjacency_matrix(adjacency_matrix, node_features_matrix)
    end
    [adjacency_matrix, node_features_matrix]
  end

  def print_adjacency_matrix(adjacency_matrix, node_features_matrix)
    i = 0
    n = adjacency_matrix.length
    (0..n - 1).each do |row|
      connectivity_index = 0 # need to set back to zero for each new row
      edge_count = 0
      neighbours = ''
      print "\n["
      (0..n - 1).each do |column|
        print adjacency_matrix[row][column]
        if (column < n - 1)
          print ', '
        else
          print ']'
        end
        connectivity_index = connectivity_index + adjacency_matrix[row][column] * (column + 1) # column index number * matrix element (0 or 1)
        if (adjacency_matrix[row][column] == 1)
          neighbours = neighbours + column.to_s + ","
          edge_count += 1
        end
      end
      neighbours.chop!
      print " [#{node_features_matrix[i]}] [#{edge_count}] [#{connectivity_index}] [#{neighbours}]\n"
      i = i + 1
    end
  end

  def write_ninchi_string(adjacency_matrix, node_features_matrix, periodic_table_elements)
    sum_formula = write_sum_formula_string(node_features_matrix, periodic_table_elements)
    serialized_molecule = serialize_molecule(adjacency_matrix)
    "nInChI=1S/#{sum_formula}/c#{serialized_molecule}"
  end

  def write_dot_file(adjacency_matrix, node_features_matrix, filename, periodic_table_elements, periodic_table_colors)
    filename = File.basename(filename, '.mol')
    dotfile = "graph #{filename}\n{\n  bgcolor=grey\n"
    n = adjacency_matrix.length
    (0..n - 1).each do |i|
      symbol = periodic_table_elements[node_features_matrix[i] - 1]
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

  def serialize_molecule(adjacency_matrix)
    #
    # inchi_string[0] -> "2-tuple format" (0-2)(1-2)
    # inchi_string[1] -> "n-tuple format" (2:0,1)
    # inchi_string[2] -> binary format with each row of adjacency matrix elements separated by minus sign
    # inchi_string[3] -> binary number generated from inchi_string[2] by removal of all "-" and addition of leading "1"
    # inchi_string[4] -> decimal format of binary inchi_string[3]
    # inchi_string[5] -> hexadecimal format of binary inchi_string[3]
    # inchi_string[6] -> base32 encoding of binary inchi_string[3]
    #
    inchi_string = Array.new(7, '')
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
    unique_elements = node_features_matrix.uniq
    initial_counts = Array.new(unique_elements.length, 0)
    element_counts = unique_elements.zip(initial_counts).to_h
    node_features_matrix.each { |atom| element_counts[atom] += 1 }
    element_counts.transform_keys! { |k| periodic_table_elements[k - 1] } # change atomic mass to element symbol; k - 1; convert from mass back to index
  end
end
