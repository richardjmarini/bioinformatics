#!/usr/bin/env python

from sys import stdin, stdout
from pprint import pprint
from nucleotide import Nucleotide

class Cell:

    def __init__(self, value, previous= None):

        self.value= value
        self.previous= previous

    def __str__(self):

        return self.value
    

class Alignment:
 
    def __init__(self, streams= 2):

         self.streams= streams
         self.buffer= []

         for stream in range(0, streams):
             self.buffer.append([])

    def feed(self, dna, stream= 0): 

         dna= map(Nucleotide, list(dna))

         self.buffer[stream]+= dna
   
    @property
    def sequence1(self):

        return ''.join([nucleotide.code for nucleotide in self.buffer[0]])

    @property
    def sequence2(self):

        return ''.join([nucleotide.code for nucleotide in self.buffer[1]])

    @staticmethod
    def get_needleman_wunsch_matrix(sequence1, sequence2, gap_score= -1):

        # create blank matrix
        matrix= []
        for i in range(0, -len(sequence2) - 1, gap_score):
            new= []
            for j in range(0, -len(sequence1) - 1, gap_score):
                new.append(0)
            matrix.append(new)

        # initalize matrix w/ gap scores
        for i in range(0, len(sequence2) + 1):
            matrix[i][0]= -i
        for j in range(0, len(sequence1) + 1):
            matrix[0][j]= -j 

        for i in range(1, len(sequence2) + 1):
            for j in range(1, len(sequence1) + 1): 

                horz= matrix[i-1][j] + gap_score
                vert= matrix[i][j-1] + gap_score

                if sequence2[i-1].code == sequence1[j-1].code:
                    diag= matrix[i-1][j-1] + 1
                else:
                    diag= matrix[i-1][j-1] - 1

                matrix[i][j]= max(horz, vert, diag)

        return matrix

    def align(self):

        self.needleman_wunsch_matrix= self.get_needleman_wunsch_matrix(self.buffer[0], self.buffer[1])

        pprint(self.needleman_wunsch_matrix)


    def __str__(self):

        m= self.needleman_wunsch_matrix
        s= ''
        s+=  '        ' + '   '.join([b.code for b in self.buffer[0]]) + '\n'
        for i in range(1, len(self.buffer[1]) + 1):
            s+= [b.code for b in self.buffer[0]][i-1]
            for j in range(1, len(self.buffer[0]) + 1):
                   s+= '%4i' % (m[i][j])
            s+= '\n'

        return s
                



if __name__ == '__main__':

    alignment= Alignment(streams= 2)
    alignment.feed('CGTGAATTCAT', 0)
    alignment.feed('GACTTAC', 1)
    alignment.align()

    #print alignment
    
