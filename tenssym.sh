#!/bin/bash
#
# Bash script for running the Python program for computation of nonzero
# and independent tensor elements for any point-symmetry group, tensor
# type (polar or axial) and tensor rank.
#
# Copyright (C) 2023, Fredrik Jonsson, under GPL 3.0. See enclosed LICENSE.
#
PYTHON="python3"
catdir="atlas"      # The name of the catalogue to be generated

#
# Define the tensor ranks for which the atlas is to be generated. The min/max
# construction is something which in the future should be altered to instead
# get the values directly from the tensorRanks variable.
#
tensorRanks=( 2 3 4 5 )
minTensorRank=2
maxTensorRank=5

#
# Define the types of tensors which for each rank should be considered.
# (There only exists two types, polar or axial.)
#
tensorType=(
    "polar"   # Polar ("true tensor")
    "axial"   # Axial ("pseudotensor")
)

#
# Define the sets of point-symmetry groups which should be used in the
# generation of the atlas. In practice, this should exactly match the
# complete list of the crystallographic point-symmetry groups as supported
# by the TENSSYM Python program.
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

plain_test() {
    $PYTHON tenssym.py --symmetry trigonal_3 --type polar --rank 4
}

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
    echo "The following compilation of nonzero elements of polar and axial"\
         "tensor elements was generated using the TensSym package"\
	 "(https://github.com/hp35/tenssym/) and covers tensors of rank"\
	 "$minTensorRank up to rank $maxTensorRank, for the following"\
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
    echo "1.   Susceptibility Tensors for Nonlinear Optics, S. V. Popov, Yu."\
         "P. Svirko and N. Zheludev (Institute of Physics Publishing,"\
         "Bristol, 1995). ISBN 0-7503-0253-4." >> $rfile
    echo "2. The Elements of Nonlinear Optics, P. N. Butcher and D. Cotter"\
         "(Cambridge University Press, Bristol, 1993). ISBN 0-521-42424-0."\
         >> $rfile
    echo "3. The Nonlinear Optics of Magneto-Optic Media (doctoral thesis),"\
         "Fredrik Jonsson (The Royal Institute of Technology, 2000),"\
         "ISBN 91-7170-575-9; Available via the National Library of Sweden"\
         "(Kungl. Biblioteket) at"\
         "http://urn.kb.se/resolve?urn=urn:nbn:se:kth:diva-2967" >> $rfile
    echo "4. Patrick Steglich and Achim Kehrein, <em>Light propagation in"\
         "anisotropic materials and electro-optical effects: tutorial on the"\
         "use of eigenvalue problems, tensors, and symmetries</em>, J. Opt."\
         "Soc. Am. B, Vol. 41, No. 9, pp. 2191--2210 (2024)." >> $rfile
    echo "5. Herbert Goldstein, <em>Classical Mechanics</em>, 2nd Edn."\
         "(Addison-Wesley, 1980), Ch. 4.4, pp. 146--147." >> $rfile
    echo "" >> $rfile
    echo "## Copyright" >> $rfile
    echo "Copyright (C) 2023, Fredrik Jonsson, under GPL 3.0."\
	 "See enclosed LICENSE." >> $rfile
    echo "" >> $rfile
    echo "## Location of master source code" >> $rfile
    echo "The source and documentation can be found at"\
         "https://github.com/hp35/tenssym " >> $rfile
    echo "This README.md was auto-generated by the Bash shell script"\
         "https://github.com/hp35/tenssym/blob/main/tenssym.sh"\
         "on `date`." >> $rfile
}

#
# Function for auto-generation of the various-rank's README.md of the
# tensor atlas.
#
generate_atlas_rank_readmes() {
    for r in "${tensorRanks[@]}"; do
        directory="rank-$r"
        rfile="$catdir/rank-$r/README.md"
        echo "# Nonzero elements of rank-$r polar and axial tensors" > $rfile
        echo "" >> $rfile
        echo "The following compilation of nonzero elements of polar (true)"\
	     "and axial (pseudo) rank-$r tensor elements was generated using"\
	     "the TensSym package (https://github.com/hp35/tenssym/)." >> $rfile
        for t in "${tensorType[@]}"; do
            subdirectory="$t"
            echo "- Rank-$r $t tensor: [$subdirectory]($subdirectory)" >> $rfile
        done
        echo "" >> $rfile
        echo "See https://github.com/hp35/tenssym for details and code"\
	     "for the implementation of the TensSym package, and"\
	     "https://github.com/hp35/tenssym/tree/main/atlas/ for the"\
	     "complete atlas of nonzero and independent elements of other"\
	     "ranks." >> $rfile
        echo "" >> $rfile
        echo "## Copyright" >> $rfile
        echo "Copyright (C) 2023, Fredrik Jonsson, under GPL 3.0."\
	     "See enclosed LICENSE." >> $rfile
        echo "" >> $rfile
        echo "## Location of master source code" >> $rfile
        echo "The source and documentation can be found at"\
             "https://github.com/hp35/tenssym " >> $rfile
        echo "This README.md was auto-generated by the Bash shell script"\
             "https://github.com/hp35/tenssym/blob/main/tenssym.sh"\
             "on `date`." >> $rfile
    done
}

#
# Function for auto-generation of the various-rank's README.md for axial
# and polar types of the tensor atlas. These README:s are located furthest
# out in the branches of the catalogue structure of the atlas of tensor
# elements, listing the point-symmetry groups given a rank and type (polar
# or axial).
#
generate_atlas_type_readmes() {
    for r in "${tensorRanks[@]}"; do
        for t in "${tensorType[@]}"; do
            rfile="$catdir/rank-$r/$t/README.md"
            echo "- Rank-$r $t tensor: [$subdirectory]($subdirectory)" >> $rfile
            echo "# Nonzero elements of rank-$r $t tensors" > $rfile
            echo "" >> $rfile
            echo "The following compilation of nonzero elements of polar"\
		 "(true) and axial (pseudo) rank-$r tensor elements was "\
		 "generated using the TensSym package"\
		 "(https://github.com/hp35/tenssym/)." >> $rfile
            echo "" >> $rfile
            echo "The summaries of non-zero and independent elements of the"\
		 "$t tensors of rank $r are available for the following"\
		 "crystallographic point-symmetry groups:" >> $rfile
            for g in "${pointSymmetryGroups[@]}"; do
		directory="$catdir/rank-$r/$t"
                summaryfile="tensor-$g-rank-$r-$t.txt"
                md5value=($(md5sum $directory/$summaryfile))
                echo "- $g: [$summaryfile]($summaryfile)"\
		     "[MD5: $md5value]" >> $rfile
            done
            echo "" >> $rfile
	    echo "See https://github.com/hp35/tenssym for details and code"\
	         "for the implementation of the TensSym package, and"\
	         "https://github.com/hp35/tenssym/tree/main/atlas/ for the"\
	         "complete atlas of nonzero and independent elements of other"\
	         "ranks." >> $rfile
            echo "" >> $rfile
            echo "## Copyright" >> $rfile
            echo "Copyright (C) 2023, Fredrik Jonsson, under GPL 3.0."\
	         "See enclosed LICENSE." >> $rfile
            echo "" >> $rfile
            echo "## Location of master source code" >> $rfile
            echo "The source and documentation can be found at"\
                 "https://github.com/hp35/tenssym " >> $rfile
            echo "This README.md was auto-generated by the Bash shell script"\
                 "https://github.com/hp35/tenssym/blob/main/tenssym.sh"\
                 "on `date`." >> $rfile
        done
    done
}

#
# Run the functions in order. Simple as that.
#
plain_test
generate_atlas
generate_atlas_root_readme
generate_atlas_rank_readmes
generate_atlas_type_readmes
