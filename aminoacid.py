#!/usr/bin/env python

class AminoAcid:

    _table= {
        'A':('Ala', 'Alanine'),
        'C':('Cys', 'Cysteine'),
        'D':('Asp', 'Aspartic Acide'),
        'E':('Glu', 'Glutamic Acid'),
        'F':('Phe', 'Phenylalanine'),
        'G':('Gly', 'Glycine'),
        'H':('His', 'Histidine'),
        'I':('Ile', 'Isoleucine'),
        'K':('Lys', 'Lysine'),
        'L':('Leu', 'Leucine'),
        'M':('Met', 'Methionine'),
        'N':('Asn', 'Asparagine'),
        'P':('Pro', 'Proline'),
        'Q':('Gln', 'Glutamine'),
        'R':('Arg', 'Arginine'),
        'S':('Ser', 'Serine'),
        'T':('Thr', 'Threonine'),
        'V':('Val', 'Valine'),
        'W':('Trp', 'Tryptophan'),
        'Y':('Tyr', 'Tyrosine')
    }

    def __init__(self, code= None):
        self.code= code
        (self.label, self.name)= self._table.get(code.upper())
        self.next= None
        self.previous= None

    def __str__(self):
    
         return self.code
