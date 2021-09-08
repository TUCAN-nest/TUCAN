require './inchi'
require 'optparse'
require './periodic_table'

# run as `ruby cli_app.rb --molfile=<path/to/molefile>`

options = {}
OptionParser.new do |opt|
  opt.on('--molfile MOLFILE') { |o| options[:molfile] = o }
end.parse!

puts 'A new International Chemical Identifier (nInChI)'
puts 'CC BY-SA | Ulrich Schatzschneider | Universität Würzburg | NFDI4Chem | v1.4 | 06/2021'

periodic_table_elements = PeriodicTable::Elements
periodic_table_colors = PeriodicTable::ElementColor

filename = options[:molfile]
molfile_data = read_molfile(filename)
puts "\nPrinting molfile: #{filename}. First 4 lines contain header."
molfile_data.each { |line| puts line }

molecule = create_molecule_array(molfile_data, periodic_table_elements)
canonicalized_molecule = canonicalize_molecule(molecule)
puts "\n#{write_ninchi_string(canonicalized_molecule, periodic_table_elements)}"
puts "\n#{write_dot_file(canonicalized_molecule, periodic_table_elements, periodic_table_colors)}"
puts 'Output format: DOT file - to display go to https://dreampuf.github.io/GraphvizOnline/#'
