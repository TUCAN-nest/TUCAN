require "./inchi_v1.4"


# label1 = FXLabel.new(self, "File: none")
# FXLabel.new(frame1, "Enter filename:")
# enterFilename = FXTextField.new(frame1, 25)
# molfile = FXText.new(frame2)
# outputFile = FXText.new(frame3)
# label2 = FXLabel.new(self, "Sum formula: none")
# label3 = FXLabel.new(self, "Output format: none")
# loadButton = FXButton.new(frame4, "Load molfile *.mol")
# nInChIButton = FXButton.new(frame4, "Create nInChI string")
# outputFileButton = FXButton.new(frame4, "Create DOT file")


def canonicalize_molecule(molecule)
    if(molecule.length > 1)
        for i in 0..(molecule.length-1)*20
            print "=== Start Pass #",i," ===\n\n"
            canonicalized_molecule = canonicalization1(molecule)
            canonicalized_molecule = canonicalization2(canonicalized_molecule)
            print "=== End Pass #",i," ===\n\n"
        end
    end
    return canonicalized_molecule
end


def create_ninchi_string(molecule)
    sumFormulaString=calculate_sum_formula(molecule)
    ninchi_string = "nInChI=1S/#{sumFormulaString}/c#{serialization(molecule)}"
    return ninchi_string
end


puts "A new International Chemical Identifier (nInChI)"
puts "CC BY-SA | Ulrich Schatzschneider | Universität Würzburg | NFDI4Chem | v1.4 | 06/2021"

puts "Please enter a file name:"
filename = gets.chomp
file = read_molfile(filename)

molecule = create_molecule_array(file,filename)
canonicalized_molecule = canonicalize_molecule(molecule)
ninchi_string = create_ninchi_string(canonicalized_molecule)
dot_file_string = create_dot_file(canonicalized_molecule)
puts("Output format: DOT file - to display go to https://dreampuf.github.io/GraphvizOnline/#")
