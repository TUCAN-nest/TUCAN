require './inchi_v1.4'
require 'optparse'
require './periodic_table'

# run as `ruby cli_app.rb --molfile=<path/to/molefile>`

options = {}
OptionParser.new do |opt|
  opt.on('--molfile MOLFILE') { |o| options[:molfile] = o }
end.parse!

# label2 = FXLabel.new(self, "Sum formula: none")
# label3 = FXLabel.new(self, "Output format: none")
# loadButton = FXButton.new(frame4, "Load molfile *.mol")
# nInChIButton = FXButton.new(frame4, "Create nInChI string")
# outputFileButton = FXButton.new(frame4, "Create DOT file")

puts 'A new International Chemical Identifier (nInChI)'
puts 'CC BY-SA | Ulrich Schatzschneider | Universität Würzburg | NFDI4Chem | v1.4 | 06/2021'

periodic_table_elements = PeriodicTable::Elements
periodic_table_colors = PeriodicTable::ElementColor

filename = options[:molfile]
molfile_data = read_molfile(filename)

molecule = create_molecule_array(molfile_data)
canonicalized_molecule = canonicalize_molecule(molecule, periodic_table_elements)
create_ninchi_string(canonicalized_molecule, periodic_table_elements)
create_dot_file(canonicalized_molecule, periodic_table_colors)
puts 'Output format: DOT file - to display go to https://dreampuf.github.io/GraphvizOnline/#'
