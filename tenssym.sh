#!/bin/bash
#
# Bash script for running the Python program for computation of nonzero
# and independent tensor elements for any point-symmetry group, tensor
# type (polar or axial) and tensor rank.
#
# Copyright (C) 2023, Fredrik Jonsson, under GPL 3.0. See enclosed LICENSE.
#
PYTHON="python3"
catdir="atlas"        # The catalogue of tensor elements to be generated

#
# Complete listing of the crystallographic point-symmetry groups
# supported by the TENSSYM Python program.
#
ppointSymmetryGroups=(
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

pointSymmetryGroups=(
    "cubic_4bar3m"    # 43m cubic, inversion symmetry
    "cubic_23"        # 23 cubic, no inversion symmetry
)

tensorType=(
    "polar"   # Polar ("true tensor")
    "axial"   # Axial ("pseudotensor")
)

ptensorRanks=( 2 3 4 5 )
tensorRanks=( 2 3 )
minTensorRank=2
maxTensorRank=5

$PYTHON tenssym.py --symmetry trigonal_3 --type polar --rank 4

#
# Function for the automatic generation of the complete atlas of non-zero
# and independent tensor elements of the ranks and point-symmetry groups
# listed by the tensorRanks and pointSymmetryGroups variables defined above.
#
generate_atlas() {
    mkdir -p "$catdir"
    for r in "${tensorRanks[@]}"; do
        mkdir -p "$catdir/rank-$r"
        for t in "${tensorType[@]}"; do
            mkdir -p "$catdir/rank-$r/$t"
        done
    done
    for r in "${tensorRanks[@]}"; do
        for t in "${tensorType[@]}"; do
            for g in "${pointSymmetryGroups[@]}"; do
                directory="$catdir/rank-$r/$t"
                file="$directory/tensor-$g-rank-$r-$t.txt"
                rm -Rf $file
                echo "Summary $file"|tee $file
                echo "Summary generated $(date) by TensSym for"|tee $file
                echo "point-symmetry group $g, $t tensor of rank $r."|tee -a $file
                echo "See https://github.com/hp35/tenssym for details."|tee -a $file
                $PYTHON tenssym.py --symmetry $g --type $t --rank $r |tee -a $file
            done
        done
    done
}

#
# Function for auto-generation of the root README.md of the tensor atlas.
#
generate_atlas_root_readme() {
    rfile="$catdir/README.md"
    echo "# Atlas of nonzero elements of polar and axial tensors" > $rfile
    echo "" >> $rfile
    echo "The following compilation of nonzero elements of polar and axial "\
         "tensor elements was generated using the TensSym package "\
	 "(https://github.com/hp35/tenssym/) and covers tensors of rank "\
	 "$minTensorRank up to rank $maxTensorRank, for the following "\
	 "point-symmetry groups:" >> $rfile
    for g in "${pointSymmetryGroups[@]}"; do
        echo "- $g" >> $rfile
    done
    echo "" >> $rfile
    echo "See https://github.com/hp35/tenssym for details." >> $rfile
    echo "" >> $rfile
    echo "## Catalogue of nonzero and independent tensor elements" >> $rfile
    for r in "${tensorRanks[@]}"; do
        directory="rank-$r"
        echo "- Rank $r tensors: [$directory]($directory)" >> $rfile
        for t in "${tensorType[@]}"; do
            subdirectory="rank-$r/$t"
            echo "    - $t tensor: [$subdirectory]($subdirectory)" >> $rfile
        done
    done
    echo "" >> $rfile
    echo "## References" >> $rfile
    echo "1.   Susceptibility Tensors for Nonlinear Optics, S. V. Popov, Yu. "\
         "P. Svirko and N. Zheludev (Institute of Physics Publishing, "\
         "Bristol, 1995). ISBN 0-7503-0253-4." >> $rfile
    echo "2. The Elements of Nonlinear Optics, P. N. Butcher and D. Cotter"\
         "(Cambridge University Press, Bristol, 1993). ISBN 0-521-42424-0."\
         >> $rfile
    echo "3. The Nonlinear Optics of Magneto-Optic Media (doctoral thesis),"\
         "Fredrik Jonsson (The Royal Institute of Technology, 2000),"\
         "ISBN 91-7170-575-9; Available via the National Library of Sweden "\
         "(Kungl. Biblioteket) at "\
         "http://urn.kb.se/resolve?urn=urn:nbn:se:kth:diva-2967" >> $rfile
    echo "4. Patrick Steglich and Achim Kehrein, <em>Light propagation in "\
         "anisotropic materials and electro-optical effects: tutorial on the "\
         "use of eigenvalue problems, tensors, and symmetries</em>, J. Opt. "\
         "Soc. Am. B, Vol. 41, No. 9, pp. 2191--2210 (2024)." >> $rfile
    echo "5. Herbert Goldstein, <em>Classical Mechanics</em>, 2nd Edn. "\
         "(Addison-Wesley, 1980), Ch. 4.4, pp. 146--147." >> $rfile
    echo "" >> $rfile
    echo "## Copyright" >> $rfile
    echo "Copyright (C) 2023, Fredrik Jonsson, under GPL 3.0. "\
	 "See enclosed LICENSE." >> $rfile
    echo "" >> $rfile
    echo "## Location of master source code" >> $rfile
    echo "The source and documentation can be found at "\
	 "https://github.com/hp35/tenssym" >> $rfile
}

#
# Run the functions in order. Simple as that.
#
generate_atlas
generate_atlas_root_readme
