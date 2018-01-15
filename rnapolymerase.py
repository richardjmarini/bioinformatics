#!/usr/bin/env python

from sys import stdin, stdout

class RNAPolymarse:

    _table= {
        'G': 'C',
        'C': 'G',
        'A': 'T',
        'T': 'A'
    }
 
    def __init__(self, dna= None):

         self.buffer= []
         if dna:
             self.feed(dna)

    def feed(self, dna): 

         self.buffer+= list(dna)

    @property
    def mRNA(self):

        while len(self.buffer):
            nucleotide= self.buffer.pop(0) 
            template= self._table.get(nucleotide)
            yield self._table.get(template).replace('T', 'U')

if __name__ == '__main__':

    rnapolymarse= RNAPolymarse()
    while True:
        nucleotides= stdin.read()
        if not nucleotides:
            break
        nucleotides= nucleotides[:-1]
        rnapolymarse.feed(nucleotides)

    print >> stdout, ''.join(rnapolymarse.mRNA)
