#!/usr/bin/env python

from argparse import ArgumentParser
from sys import stdin, stdout

from aminoacid import AminoAcid

class Enzyme:

    _table= {
        "BamHI": {
            "name": "Bacillus amyloliquefaciens H",
            "restriction_sites": {
                '5': "C|GATCC",
                '3':  "CCTAG|G"
            }
        },
        "BAlI": {
            "name": "Brevibacterium albidum",
            "restriction_sites": {
                '5': "TGG|CCA", 
                '3': "ACC|GGT"
            }
        },
        "EcoRI": {
            "name": "Escherichia coli RY13",
            "restriction_sites": {
                '5': "G|AATTC",
                '3': "CTTAA|G"
            }
        },
        "HaeII": {
            "name": "Haemophilus aegyptuis",
            "restriction_sites": {
                '5': "[AG]GGC|[TC]", 
                'c': "[TC]|CGGC[AG]"
            }
        },
        "HaeIII": {
            "name": "Haemophilus aegyptuis", 
            "restriction_sites": {
                '5': "GG|CC",
                '3': "CC|GG"
            }
        },
        "HindII": {
            "name": "Haemophilus influenzae Rd",
            "restriction_sites": {
                '5': "GT[TC]|[AG]AC", 
                '3': "CA[AG]|[TC]TG"
            }
        },
        "HindIII": {
            "name": "Haemophilus influenzae Rd",
            "restriction_sites": {
                '5': "A|AGCTT", 
                '3': "TTCGA|A"
            }
        },
        "HpaI": {
            "name": "Haemophilius parainfluenzae",
            "restriction_sites": {
                '5': "GTT|AAC", 
                '3': "GAA|TTC"
            }
        },
        "HpaII": {
            "name": "Haemophilius parainfluenzae",
            "restriction_sites": {
                '5': "C|CGG",
                '3': "GGC|C"
            }
        },
        "psTI": {
            "name": "Providencia stuartii 164",
            "restriction_sites": {
                '5': "CTGCA|A",
                '3': "C|ACGTC"
            }
        },
        "SalI": {
            "name": "Streptomyces albus G",
            "restriction_sites": {
                '5': "G|TCGAC",
                '3': "CAGCT|G"
            }
        }
    }

    def __init__(self, code):
        self.code= code
        self.__dict__.update(**self._table.get(self.code))


if __name__ == '__main__':

    parser= ArgumentParser(description= 'restricts nucleotide sequence(s)')
    parser.add_argument('--sequences', nargs= '+', help= 'list of nucleotide sequences')
    parser.add_argument('--enzymes', nargs= '+', help= 'list of enzymes to apply to nucleotide sequence(s)')

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

    # TODO: 
    # build restriction map class
    # create activate() method
