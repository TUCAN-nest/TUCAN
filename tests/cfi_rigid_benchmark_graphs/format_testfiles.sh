#!/bin/bash

# Extract all files from https://www.lics.rwth-aachen.de/go/id/rtok/.
tar -xvf cfi-rigid-z2-tar.gz --wildcards '*-1' --strip-components 1
tar -xvf cfi-rigid-r2-tar.gz --wildcards '*-1' --strip-components 1
tar -xvf cfi-rigid-s2-tar.gz --wildcards '*-1' --strip-components 1
tar -xvf cfi-rigid-t2-tar.gz --wildcards '*-1' --strip-components 1
tar -xvf cfi-rigid-z3-tar.gz --wildcards '*-1' --strip-components 1
tar -xvf cfi-rigid-d3-tar.gz --wildcards '*-1' --strip-components 1

# Remove all files that contain graphs with more than 1000 nodes. The number of
# nodes is indicated in position 15 through 18 in the filename.
for f in ./*-1; do
    if [ ${f:15:4} -gt 1000 ]; then
        echo ${f:15:4}
        rm $f
    fi
done

# Append suffix.
for f in ./*-1; do
    mv $f "${f}.col"
done
