#!/bin/bash
#
# Bash script for running the Python program for computation of nonzero
# and independent tensor elements for any point-symmetry group, tensor
# type (polar or axial) and tensor rank.
#
# Copyright (C) 2023, Fredrik Jonsson. Non-commercial copying welcome.
#
PYTHON="python3"
catdir="catalogue"    # The catalogue of tensor elements to be generated

#
# Complete listing of the crystallographic point-symmetry groups
# supported by the TENSSYM Python program.
#
pointSymmetryGroups=(
    "cubic_4bar3m"    # 43m cubic, inversion symmetry
    "cubic_23"        # 23 cubic, no inversion symmetry
    "trigonal_3"      # 3 trigonal, no inversion symmetry
    "trigonal_32_2x"  # 32 trigonal, no inversion symmetry, 2-fold x-symmetry
    "trigonal_32_2y"  # 32 trigonal, no inversion symmetry, 2-fold y-symmetry
    "triclinic_1bar"  # 1 triclinic, inversion symmetry for all axes
    "triclinic_1"     # 1 triclinic, no inversion symmetry
    "tetragonal_422"  # 422 tetragonal, no inversion symmetry
    "cubic_m3m"       # m3m cubic, inversion symmetry
    "cubic_432"       # 432 cubic, no inversion symmetry
    "cubic_m3"        # m3 cubic, inversion symmetry
    "hexagonal_6"     # 6 hexagonal, no inversion symmetry
    "hexagonal_6bar"  # 6 hexagonal, inversion symmetry <= GER FEL ANTAL TERMER!
)

tensorType=(
    "polar"   # Polar ("true tensor")
    "axial"   # Axial ("pseudotensor")
)

tensorRank=( 2 3 4 5 )

$PYTHON tenssym.py --symmetry trigonal_3 --type polar --rank 4

mkdir "$catdir"
for k in "${tensorRank[@]}"; do
    mkdir "$catdir/rank-$k"
    for j in "${tensorType[@]}"; do
        mkdir "$catdir/rank-$k/$j"
    done
done
for i in "${pointSymmetryGroups[@]}"; do
    for j in "${tensorType[@]}"; do
        for k in "${tensorRank[@]}"; do
	    directory="$catdir/rank-$k/$j"
	    file="$directory/tensor-$i-rank-$k-$j.txt"
	    rm -Rf $file
	    echo "Summary $file" |tee $file
	    echo "Summary generated $(date) for" |tee $file
            echo "point-symmetry group $i, $j tensor of rank $k." |tee -a $file
	    $PYTHON tenssym.py --symmetry $i --type $j --rank $k |tee -a $file
        done
    done
done
