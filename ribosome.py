#!/usr/bin/env python

from argparse import ArgumentParser
from sys import stdin, stdout

from aminoacid import AminoAcid

class Ribosome:

    _table= {
        'G': {
            'G': {
                'G': AminoAcid('G'),
                'A': AminoAcid('G'),
                'C': AminoAcid('G'),
                'U': AminoAcid('G'),
            },
            'A': {
                'G': AminoAcid('E'),
                'A': AminoAcid('E'),
                'C': AminoAcid('D'),
                'U': AminoAcid('D'),
            },
            'C': {
                'G': AminoAcid('A'),
                'A': AminoAcid('A'),
                'C': AminoAcid('A'),
                'U': AminoAcid('A'),
            },
            'U': {
                'G': AminoAcid('V'),
                'A': AminoAcid('V'),
                'C': AminoAcid('V'),
                'U': AminoAcid('V'),
            }
        },
        'A': {
            'G': {
                'G': AminoAcid('R'),
                'A': AminoAcid('R'),
                'C': AminoAcid('S'),
                'U': AminoAcid('S'),
            },
            'A': {
                'G': AminoAcid('K'),
                'A': AminoAcid('K'),
                'C': AminoAcid('N'),
                'U': AminoAcid('N'),
            },
            'C': {
                'G': AminoAcid('T'),
                'A': AminoAcid('T'),
                'C': AminoAcid('T'),
                'U': AminoAcid('T'),
            },
            'U': {
                'G': AminoAcid('M'),
                'A': AminoAcid('I'),
                'C': AminoAcid('I'),
                'U': AminoAcid('I'),
            }
        },
        'C': {
            'G': {
                'G': AminoAcid('R'),
                'A': AminoAcid('R'),
                'C': AminoAcid('R'),
                'U': AminoAcid('R'),
            },
            'A': {
                'G': AminoAcid('Q'),
                'A': AminoAcid('Q'),
                'C': AminoAcid('H'),
                'U': AminoAcid('H'),
            },
            'C': {
                'G': AminoAcid('P'),
                'A': AminoAcid('P'),
                'C': AminoAcid('P'),
                'U': AminoAcid('P'),
            },
            'U': {
                'G': AminoAcid('L'),
                'A': AminoAcid('L'),
                'C': AminoAcid('L'),
                'U': AminoAcid('L'),
            }
        },
        'U': {
            'G': {
                'G': AminoAcid('W'),
                'A': 'STOP',
                'C': AminoAcid('C'),
                'U': AminoAcid('C'),
            },
            'A': {
                'G': 'STOP',
                'A': 'STOP',
                'C': AminoAcid('Y'),
                'U': AminoAcid('Y'),
            },
            'C': {
                'G': AminoAcid('S'),
                'A': AminoAcid('S'),
                'C': AminoAcid('S'),
                'U': AminoAcid('S'),
            },
            'U': {
                'G': AminoAcid('L'),
                'A': AminoAcid('L'),
                'C': AminoAcid('F'),
                'U': AminoAcid('F'),
            }
        }
    }


    def __init__(self, mRNA= None):

         self.buffer= []
         if mRNA:
             self.feed(mRNA)

    def feed(self, mRNA): 

         self.buffer+= list(mRNA)

    @property
    def peptides(self):

        node= self._table
        while len(self.buffer):

            nucleotide= self.buffer.pop(0)
            node= node.get(nucleotide)
            if isinstance(node, AminoAcid):
                yield node
                node= self._table

            elif node == 'STOP': 
                yield 
                node= self._table

    @property
    def peptide_chains(self):
        chain= []
        for peptide in self.peptides:
            if peptide: 
                chain.append(peptide)
            else:
                yield tuple(chain)
                chain= []
        yield chain


if __name__ == '__main__':

    parser= ArgumentParser(description ='Aligns nucleotide sequences')
    parser.add_argument('--sequences', nargs= '+', help= 'list of nucleotide sequences')

    arguments= parser.parse_args()
    if arguments.sequences is None:
        print "missing sequences argument"
        parser.print_help()
        exit(-1)

    ribosome= Ribosome()

    if '-' in arguments.sequences:

        while True:
             nucleotides= stdin.read()
             if not nucleotides:
                 break
             nucleotides= nucleotides[:-1]
             ribosome.feed(nucleotides)
    else:
        map(ribosome.feed, arguments.sequences)

    for chain in ribosome.peptide_chains:
        print >> stdout, ''.join(map(str, chain))
