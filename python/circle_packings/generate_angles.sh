#! /bin/bash

let a=$1   # number of vertices to enumerate from
let a=$a+2 # convert to number of faces
let c=$2   # number of graphs to skip in order of the output from Plantri (given a)
let sum=0

while test 0 -eq 0
do
    rm -f ./graphs
    
    touch ./graphs
    
    plantri -qc4m3d $a >> ./graphs 2>tmp
    
    b=`grep quartic tmp | (read a junk; echo $a)`
    
    let d=$a-2
    
    if ./generate_angles.py $d $c $b
        then
        :
        else
        break
    fi
    
    let c=0
    
    let sum=$sum+$b
    
    let a=$a+1
done

rm -f tmp

echo "Enumerated ${sum} angle structures (with volumes)."

# -q for planar quadrangulations
# -c4 for 3-connected and no 6-cyclically-edge-connected dual
# -m3 for minimum degree 3 (no bigons in dual - possibly already guaranteed and hence not necessary)
# -o for enumeration up to orientation preserving isomorphism
# -a for ASCII output
# -d for duals (lowest precedence - applied directly before outputting)
