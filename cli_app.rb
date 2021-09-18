require './inchi'
require 'optparse'

# run as `ruby cli_app.rb --molfile=<path/to/molefile>`
# run with `--permute-input` flag in order to randomly permute the atom indices

class CommandLineInterface
  include Inchi

  def initialize
    options = {}
    OptionParser.new do |opt|
      opt.on('--molfile MOLFILE') { |o| options[:molfile] = o }
      opt.on('--permute-input') { |o| options[:permute] = o }
    end.parse!
    @filename = options[:molfile]
    @permute = options[:permute] || false
  end

  def run
    puts "\nA new International Chemical Identifier (nInChI)"
    puts 'CC BY-SA | Ulrich Schatzschneider | Universität Würzburg | NFDI4Chem | v1.4 | 06/2021'

    molfile_data = read_molfile(@filename)
    puts "\nPrinting molfile: #{@filename}. First 4 lines contain header."
    molfile_data.each { |line| puts line }

    molecule = create_molecule_array(molfile_data)
    molecule = update_molecule_indices(molecule, random_indices: true) if @permute
    canonicalized_molecule = canonicalize_molecule(molecule, @filename)
    puts "\nnInChI string for #{File.basename(@filename, '.mol')}:\n#{write_ninchi_string(canonicalized_molecule)}"
    puts "\nDOT file for #{File.basename(@filename, '.mol')} (display graph at https://dreampuf.github.io/GraphvizOnline/#):\n#{write_dot_file(canonicalized_molecule, @filename)}"
    puts "#{'-' * 100}\n#{'-' * 100}"
  end
end

CommandLineInterface.new.run
