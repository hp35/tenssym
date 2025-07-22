#
# Makefile to be used for the generation of a catalogue of non-zero and
# independent tensor elements of polar and axial tensors.
#
# Copyright (C) 2023, Fredrik Jonsson, under GPL 3.0. See enclosed LICENSE.
#

#
# Re-generate the entire atlas of non-zero and independent tensor elements, as
# present in the repository at https://github.com/hp35/tenssym/tree/main/atlas
#
all:
	./tenssym.sh

rank-4-polar-43m-rotated:
	python3 tenssym.py --symmetry trigonal_3 --type polar --rank 4 \
		--phi 0.3 --theta 0.2 --psi 0.1

#
# Generate a change log of the history of the CWEB master.
#
history:
	git log -- filename tenssym.py > HISTORY.txt

#
# Basic clean-up.
#
clean:
	-rm -Rf *~
