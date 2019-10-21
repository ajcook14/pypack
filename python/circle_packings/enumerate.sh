#! /bin/bash

rm -f /python/circle_packings/graphs
rm -f ./polyhedra/*

touch ./graphs
let a=8
let n=$1+2
while test $a -le $n
do
    ../../plantri/plantri -qc4m3d $a >> ./graphs #> /dev/null 2>&1
    let a=$a+1
done

# -q for planar quadrangulations
# -c4 for 3-connected and no 6-cyclically-edge-connected dual
# -m3 for minimum degree 3 (no bigons in dual - possibly already guaranteed and hence not necessary)
# -o for enumeration up to orientation preserving isomorphism
# -a for ASCII output
# -d for duals (lowest precedence - applied directly before outputting)
