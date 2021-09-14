require './inchi'
require 'optparse'
require './periodic_table'

# run as `ruby cli_app.rb --molfile=<path/to/molefile>`

class CommandLineInterface
  include Inchi
  include PeriodicTable

  def initialize
    options = {}
    OptionParser.new do |opt|
      opt.on('--molfile MOLFILE') { |o| options[:molfile] = o }
    end.parse!
    @filename = options[:molfile]
  end

  def run
    puts "\nA new International Chemical Identifier (nInChI)"
    puts 'CC BY-SA | Ulrich Schatzschneider | Universität Würzburg | NFDI4Chem | v1.4 | 06/2021'

    molfile_data = read_molfile(@filename)
    puts "\nPrinting molfile: #{@filename}. First 4 lines contain header."
    molfile_data.each { |line| puts line }

    molecule = create_molecule_array(molfile_data, PeriodicTable::ELEMENTS)
    canonicalized_molecule = canonicalize_molecule(molecule)
    puts "\n#{write_ninchi_string(canonicalized_molecule, PeriodicTable::ELEMENTS)}"
    puts "\n#{write_dot_file(canonicalized_molecule, PeriodicTable::ELEMENTS, PeriodicTable::ELEMENT_COLORS)}"
    puts 'Output format: DOT file - to display go to https://dreampuf.github.io/GraphvizOnline/#'
  end
end

CommandLineInterface.new.run
