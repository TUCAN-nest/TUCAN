#
# (c) CC BY-SA | Jan Brammer, RWTH Aachen and Ulrich Schatzschneider, Universität Würzburg | NFDI4Chem | v3.0.3 | 05.12.2021
#

require 'optparse'
require './periodic_table'
require './inchi'

#
# main program loop
#
# run as `ruby cli_app.rb --molfile=<path/to/molefile>`
# run with `--permute-input` flag in order to randomly permute the atom indices
# 

class CommandLineInterface
  include PeriodicTable

  def initialize
    options = {}
    begin
      OptionParser.new do |opts|
        opts.on('--molfile MOLFILE') { |o| options[:molfile] = o }
        opts.on('--permute-input') { |o| options[:permute_input] = o }
        opts.on('--print-dotfile') { |o| options[:print_dotfile] = o }
        opts.on('--print-molfile') { |o| options[:print_molfile] = o }
      end.parse!
    rescue OptionParser::InvalidOption => e
      abort("#{e}")
    end
    @filename = options[:molfile]
    abort('Please provide a filename.') if @filename.nil?
    abort("File `#{@filename}` doesn't exist.") unless File.exist?(@filename)
    @permute_input = options[:permute_input] || false
    @print_dotfile = options[:print_dotfile] || false
    @print_molfile = options[:print_molfile] || false
  end

  def run
    puts "#{'-' * 100}\n"
    puts "\nA new International Chemical Identifier (nInChI) v3.0.3\n"
    puts "\nJan Brammer (RWTH Aachen) and Ulrich Schatzschneider (Universität Würzburg) within NFDI4Chem\n"
    puts "\nCC BY-SA 12/2021\n"
    puts "\n#{'-' * 100}\n"
    molfile_data = read_molfile(@filename)
    puts "\nPrinting molfile: #{@filename} First 6 lines contain header."
    puts "\n#{'-' * 75}\n"
    puts molfile_data
    puts "\n#{'-' * 75}\n"
    adjacency_matrix, node_features_matrix, distance_matrix, molfile_header, molfile_footer = initialize_matrix(
      molfile_data, PeriodicTable::ELEMENTS
    )
    adjacency_matrix, node_features_matrix, distance_matrix = sort_adjacency_matrix(adjacency_matrix,
                                                                                    node_features_matrix, distance_matrix)
    puts "\nFINAL STAGE\n"
    puts "\nPrinting: [adjacency matrix] [distance matrix] [node features matrix] {neighbours}\n\n"
    print_adjacency_matrix(adjacency_matrix, node_features_matrix, distance_matrix)
    puts "\n#{write_ninchi_string(adjacency_matrix, node_features_matrix, PeriodicTable::ELEMENTS)}"
    puts "\n#{write_dot_file(adjacency_matrix, node_features_matrix, @filename, PeriodicTable::ELEMENTS,
                             PeriodicTable::ELEMENT_COLORS)}" if @print_dotfile
    puts "\n#{write_molfile(adjacency_matrix, node_features_matrix, molfile_header, molfile_footer,
                            PeriodicTable::ELEMENTS)}" if @print_molfile
    puts "\nOutput format: DOT file - to display go to https://dreampuf.github.io/GraphvizOnline/#" if @print_dotfile
    puts "\n#{'-' * 100}\n"
  end
end

CommandLineInterface.new.run
