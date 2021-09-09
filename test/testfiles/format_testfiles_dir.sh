# Move molfiles, metadata, and drawings at root to a subfolder named after the molecule.

shopt -s extglob # enables use of patterns when stripping strings
shopt -s nullglob # don't execute loop in case there's no match

for moldata in  ./*.mol ./*_metadata.txt ./*_drawing.jpeg; do
  molname=${moldata%@(.mol|_metadata.txt|_drawing.jpeg)}
  if [ ! -d "$molname" ]; then
    echo "Making ${molname} directory."
    mkdir "${molname}"
  fi
  mv "$moldata" "${molname}"
done