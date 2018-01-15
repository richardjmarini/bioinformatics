RIBOSOME=./ribosome.py
RNAPOLYMERASE=./rnapolymerase.py

all: ribosome rnapolymerase protein

ribosome:
	echo 'Ribosome Test....'
	echo 'GGGGGAGGCGGU' | $(RIBOSOME)
	echo 'GAGGAAGACGAU' | $(RIBOSOME)
	echo 'GCGGCAGCCGCU' | $(RIBOSOME)
	echo 'GUGGUAGUCGUG' | $(RIBOSOME)

	echo 'AGGAGAAGCAGU' | $(RIBOSOME)
	echo 'AAGAAAAACAAU' | $(RIBOSOME)
	echo 'ACGACAACCACU' | $(RIBOSOME)
	echo 'AUGAUAAUCAUU' | $(RIBOSOME)

	echo 'CGGCGACGCCGU' | $(RIBOSOME)
	echo 'CAGCAACACCAU' | $(RIBOSOME)
	echo 'CCGCCACCCCCU' | $(RIBOSOME)
	echo 'CUGCUACUCCUU' | $(RIBOSOME)

	echo 'UGGUGAUGCUGU' | $(RIBOSOME)
	echo 'UAGUAAUACUAU' | $(RIBOSOME)
	echo 'UCGUCAUCCUCU' | $(RIBOSOME)
	echo 'UUGUUAUUCUUU' | $(RIBOSOME) 

rnapolymerase:
	echo
	echo 'RNA Polymerase Tests...'
	echo 'GACTTGAC' | $(RNAPOLYMERASE)

protein:
	echo
	echo 'DNA -> Protein Tests...'
	echo 'GACTTGAC' | $(RNAPOLYMERASE) | $(RIBOSOME)

