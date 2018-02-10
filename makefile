RIBOSOME=./ribosome.py
RNAPOLYMERASE=./rnapolymerase.py
ALIGNMENT=./alignment.py

all: ribosome rnapolymerase protein alignment

ribosome:
	echo 'Ribosome Test....'
	echo 'GGGGGAGGCGGU' | $(RIBOSOME) --sequences -
	echo 'GAGGAAGACGAU' | $(RIBOSOME) --sequences -
	echo 'GCGGCAGCCGCU' | $(RIBOSOME) --sequences -
	echo 'GUGGUAGUCGUG' | $(RIBOSOME) --sequences -

	echo 'AGGAGAAGCAGU' | $(RIBOSOME) --sequences -
	echo 'AAGAAAAACAAU' | $(RIBOSOME) --sequences -
	echo 'ACGACAACCACU' | $(RIBOSOME) --sequences -
	echo 'AUGAUAAUCAUU' | $(RIBOSOME) --sequences -

	echo 'CGGCGACGCCGU' | $(RIBOSOME) --sequences -
	echo 'CAGCAACACCAU' | $(RIBOSOME) --sequences -
	echo 'CCGCCACCCCCU' | $(RIBOSOME) --sequences -
	echo 'CUGCUACUCCUU' | $(RIBOSOME) --sequences -

	echo 'UGGUGAUGCUGU' | $(RIBOSOME) --sequences -
	echo 'UAGUAAUACUAU' | $(RIBOSOME) --sequences -
	echo 'UCGUCAUCCUCU' | $(RIBOSOME) --sequences -
	echo 'UUGUUAUUCUUU' | $(RIBOSOME) --sequences -

rnapolymerase:
	echo
	echo 'RNA Polymerase Tests...'
	echo 'GACTTGAC' | $(RNAPOLYMERASE)

protein:
	echo
	echo 'DNA -> Protein Tests...'
	echo 'GACTTGAC' | $(RNAPOLYMERASE) | $(RIBOSOME) --sequences -

alignment:
	echo
	echo 'Alignment tests....'
	echo "CGTGAATTCATGACTTAC" | $(ALIGNMENT) --sequences -

