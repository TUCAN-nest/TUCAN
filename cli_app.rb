require "./inchi_v1.4"

# label2 = FXLabel.new(self, "Sum formula: none")
# label3 = FXLabel.new(self, "Output format: none")
# loadButton = FXButton.new(frame4, "Load molfile *.mol")
# nInChIButton = FXButton.new(frame4, "Create nInChI string")
# outputFileButton = FXButton.new(frame4, "Create DOT file")


puts "A new International Chemical Identifier (nInChI)"
puts "CC BY-SA | Ulrich Schatzschneider | Universität Würzburg | NFDI4Chem | v1.4 | 06/2021"

puts "Please enter a file name:"
filename = gets.chomp
file = read_molfile(filename)

molecule = create_molecule_array(file,filename)
canonicalized_molecule = canonicalize_molecule(molecule)
ninchi_string = create_ninchi_string(canonicalized_molecule)
dot_file_string = create_dot_file(canonicalized_molecule)
puts "Output format: DOT file - to display go to https://dreampuf.github.io/GraphvizOnline/#"
