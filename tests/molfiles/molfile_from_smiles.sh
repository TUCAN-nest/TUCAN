# Create a v3000 molfile from a SMILES-string.
# First command line parameter is SMILES string (between double-quotes).
# Second command line parameter is the name of the molecule (without blanks).
smiles=$1
name=$2
mkdir "${name}"
obabel -:"$smiles" -h -omol -O "${name}/${name}.mol" -x3
