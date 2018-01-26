#!/usr/bin/env python

from sys import stdin, stdout
from pprint import pprint
from nucleotide import Nucleotide

class Alignment:
 
    def __init__(self, sequences= []):

         self.streams= len(sequences)
         self.sequences= []
         for sequence in sequences:
             self.sequences.append(map(Nucleotide, list(sequence)))

    def feed(self, dna, stream= 0): 

         dna= map(Nucleotide, list(dna))
 
         try:
             sequence= self.sequences[stream]
         except IndexError:
              for i in range(0, stream + 1):
                  self.sequences.append([])
          
         self.sequences[stream]+= dna
   
    @staticmethod
    def needleman_wunsch(sequence1, sequence2, gap_score= -1, match= 1,\
        mismatch= -1):

        # create blank matrix
        matrix= []
        for row in range(0, len(sequence2) + 1): 
            new_row= [0] * (len(sequence1) + 1)
            matrix.append(new_row)

        # initalize matrix w/ gap scores
        for row in range(0, len(sequence2) + 1):
            matrix[row][0]= gap_score * row
        for col in range(0, len(sequence1) + 1):
            matrix[0][col]= gap_score * col

        for row in range(1, len(sequence2) + 1):
            for col in range(1, len(sequence1) + 1): 

                horz= matrix[row-1][col] + gap_score
                vert= matrix[row][col-1] + gap_score

                if sequence2[row-1].code == sequence1[col-1].code:
                    diag= matrix[row-1][col-1] + match
                else:
                    diag= matrix[row-1][col-1] + mismatch

                matrix[row][col]= max(horz, vert, diag)

        return matrix

    def compute_simularity(self, algorithm= "needleman_wunsch"):

        if hasattr(self, algorithm):
            handler= getattr(self, algorithm)
            matrix= handler(self.sequences[0], self.sequences[1])
        else:
            raise Exception("Unkown Algorithm: %s" % (algorithm))

        return matrix 

    def traceback(self, matrix, row, col, alignment_s, alignment_t,\
        alignment_position= 0, gap_score= -1, match= 1, mismatch= -1):

        if row <= 0 or col <= 0:

            alignment_position= 0
 
            if row <= 0 and col > 0:
                alignment_t[alignment_position]= '-'
                alignment_s[alignment_position]= self.sequences[1][row-1].code
            elif col <= 0 and row > 0:
                alignment_t[alignment_position]= self.sequences[0][col-1].code
                alignment_s[alignment_position]= '-'
            else:
                alignment_t[alignment_position]= self.sequences[0][col-1].code
                alignment_s[alignment_position]= self.sequences[1][row-1].code

        elif row > 0 and col > 0\
            and matrix[row][col] == (matrix[row][col-1] + gap_score):

            (alignment_s, alignment_t, alignment_position)= self.traceback(\
                matrix, row, col-1, alignment_s, alignment_t,\
                alignment_position, gap_score, match, mismatch
            )

            alignment_position+= 1
            alignment_t[alignment_position]= '-'
            alignment_s[alignment_position]= self.sequences[0][col-1].code

        elif row > 0 and col > 0\
            and ( (matrix[row][col] == (matrix[row-1][col-1] + match))\
            or (matrix[row][col] == (matrix[row-1][col-1] + mismatch)) ):

            (alignment_s, alignment_t, alignment_position)= self.traceback(\
                matrix, row-1, col-1, alignment_s, alignment_t,\
                alignment_position, gap_score, match, mismatch\
             )

            alignment_position+= 1
            alignment_t[alignment_position]= self.sequences[1][row-1].code
            alignment_s[alignment_position]= self.sequences[0][col-1].code

        else: #elif row > 0 and col > 0\
            #and matrix[row][col] == (matrix[row-1][col] + gap_score):

            (alignment_s, alignment_t, alignment_position)= self.traceback(\
                matrix, row-1, col, alignment_s, alignment_t,\
                alignment_position, gap_score, match, mismatch\
            )
        
            alignment_position+= 1  
            alignment_t[alignment_position]= self.sequences[1][row-1].code
            alignment_s[alignment_position]= '-'

        return (alignment_s, alignment_t, alignment_position)


if __name__ == '__main__':

    sequence1= 'CGTGAATTCAT'
    sequence2= 'GACTTAC'

    alignment= Alignment(sequences= [sequence1, sequence2])
    #alignment.feed('CGTGA', 0)
    #alignment.feed('ATTCAT', 0)
    #alignment.feed('GACTTAC', 1)

    print 'sequences...'
    for sequence in alignment.sequences:
        for nucleotide in sequence:
            print nucleotide.code,
        print
    print

    matrix= alignment.compute_simularity("needleman_wunsch")
    print 'needleman-wunsch matrix...'
    pprint(matrix)
    print

    length_t= len(sequence1)
    length_s= len(sequence2)
    length= max(length_s, length_t)
    alignment_s= ['-'] * length
    alignment_t= ['-'] * length
   
    alignment.traceback(matrix, length_s, length_t, alignment_s, alignment_t)

    print 'alignments...'
    print alignment_s
    print alignment_t
    print
