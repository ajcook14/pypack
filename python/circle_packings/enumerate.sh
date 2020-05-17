#! /bin/bash

rm -f ./graphs
rm -f ./polyhedra/*

touch ./graphs
let a=8
let n=$1+2 # convert to number of faces
let sum=0
while test $a -le $n
do
    plantri -qc4m3d $a >> ./graphs 2>tmp
    let a=$a+1
    b=`grep quartic tmp | (read a junk; echo $a)`
    let sum=$sum+$b
done

rm tmp

echo "Enumerated ${sum} edge graphs. I am now finding the geometric structures and volumes."

python3 pack.py $sum

# -q for planar quadrangulations
# -c4 for 3-connected and no 6-cyclically-edge-connected dual
# -m3 for minimum degree 3 (no bigons in dual - possibly already guaranteed and hence not necessary)
# -o for enumeration up to orientation preserving isomorphism
# -a for ASCII output
# -d for duals (lowest precedence - applied directly before outputting)
