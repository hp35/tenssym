#
# Makefile to be used for the generation of a catalogue of non-zero and
# independent tensor elements of polar and axial tensors.
#
# Copyright (C) 2023, Fredrik Jonsson, under GPL 3.0. See enclosed LICENSE.
#
all:
	./tenssym.sh

#
# Generate a change log of the history of the CWEB master.
#
history:
	git log -- filename tenssym.py > HISTORY.txt

clean:
	-rm -Rf *~
