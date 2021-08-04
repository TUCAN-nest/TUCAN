require './inchi_v1.4'
require 'optparse'

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

filename = options[:molfile]
molfile_data = read_molfile(filename)

molecule = create_molecule_array(molfile_data)
canonicalized_molecule = canonicalize_molecule(molecule)
create_ninchi_string(canonicalized_molecule)
create_dot_file(canonicalized_molecule)
puts 'Output format: DOT file - to display go to https://dreampuf.github.io/GraphvizOnline/#'
