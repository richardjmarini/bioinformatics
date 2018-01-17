#!/usr/bin/env python

from sys import stdin, stdout
from pprint import pprint
from nucleotide import Nucleotide

class Alignment:
 
    def __init__(self, streams= 2):
         self.streams= streams
         self.buffer= []
         for stream in range(0, streams):
             self.buffer.append([])

    def feed(self, dna, stream= 0): 
         dna= map(Nucleotide, list(dna))
         self.buffer[stream]+= dna

    @staticmethod
    def get_simularity_matrix(buffer1, buffer2):

        matrix= []
        for j in range(0, len(buffer2)):
            new= []
            for i in range(0, len(buffer1)):
                new.append(0) 
            matrix.append(new)

        for i in range(0, len(buffer1)):

            for j in range(0, len(buffer2)):
                if buffer1[i].code == buffer2[j].code:
                    # nucleotide match
                    matrix[j][i]= 2
                elif buffer1[i].type == buffer2[j].type:
                    # transition:
                    # purine w/ purine or pyrimidine w/ pyrimidine
                    matrix[j][i]= 1
                else:          
                    # transversion:
                    # purine w/ pyrimidine 
                    matrix[j][i]= -1

        return matrix

    @property
    def alignment(self):

        self.simularity_matrix= Alignment.get_simularity_matrix(self.buffer[0], self.buffer[1])

        return self.simularity_matrix



if __name__ == '__main__':

    alignment= Alignment(streams= 2)
    alignment.feed('AGCG', 0)
    alignment.feed('ACGT', 1)

    alignment.alignment
    pprint(alignment.simularity_matrix, indent= 1, width= 18)
    
