#!/usr/bin/env python

from argparse import ArgumentParser
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

def activate(sequence):

    rnapolymarse= RNAPolymarse()
    rnapolymarse.feed(sequence)

    print >> stdout, ''.join(rnapolymarse.mRNA)


if __name__ == '__main__':

    parser= ArgumentParser(description= 'Converts nucleotide sequences into mRNA')
    parser.add_argument('--sequences', default= '-', nargs= '+', help= 'list of nucleotide sequences [default: -]')

    arguments= parser.parse_args()
    if arguments.sequences is None:
        print "missing sequences argument"
        parser.print_help()
        exit(-1)

    if '-' in arguments.sequences:

        while True:

            nucleotides= stdin.read()
            if not nucleotides:
                break

            nucleotides= nucleotides[:-1]
            sequences= nucleotides.split("\r")

    else:
        sequences= arguments.sequences

    map(activate, sequences)

