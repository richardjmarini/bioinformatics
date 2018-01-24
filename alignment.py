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
   
    @property
    def sequence1(self):

        return ''.join([nucleotide.code for nucleotide in self.buffer[0]])

    @property
    def sequence2(self):

        return ''.join([nucleotide.code for nucleotide in self.buffer[1]])

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
            matrix= handler(self.buffer[0], self.buffer[1])
        else:
            raise Exception("Unkown Algorithm: %s" % (algorithm))

        return matrix 

    def align(self, matrix, row, col, alignment_left, alignment_top,\
        alignment_position= 0, gap_score= -1, match= 1, mismatch= -1):

        if row <= 0 or col <= 0:

            alignment_position= 0
 
            if row <= 0 and col > 0:
                alignment_top[alignment_position]= '-'
                alignment_left[alignment_position]= self.buffer[1][row-1].code
            elif col <= 0 and row > 0:
                alignment_top[alignment_position]= self.buffer[0][col-1].code
                alignment_left[alignment_position]= '-'
            else:
                alignment_top[alignment_position]= self.buffer[0][col-1].code
                alignment_left[alignment_position]= self.buffer[1][row-1].code

        elif row > 0 and col > 0\
            and matrix[row][col] == (matrix[row][col-1] + gap_score):

            (alignment_left, alignment_top, alignment_position)= self.align(\
                matrix, row, col-1, alignment_left, alignment_top,\
                alignment_position, gap_score, match, mismatch
            )

            alignment_position+= 1
            alignment_top[alignment_position]= '-'
            alignment_left[alignment_position]= self.buffer[0][col-1].code

        elif row > 0 and col > 0\
            and ( (matrix[row][col] == (matrix[row-1][col-1] + match))\
            or (matrix[row][col] == (matrix[row-1][col-1] + mismatch)) ):

            (alignment_left, alignment_top, alignment_position)= self.align(\
                matrix, row-1, col-1, alignment_left, alignment_top,\
                alignment_position, gap_score, match, mismatch\
             )

            alignment_position+= 1
            alignment_top[alignment_position]= self.buffer[1][row-1].code
            alignment_left[alignment_position]= self.buffer[0][col-1].code

        else: #elif row > 0 and col > 0\
            #and matrix[row][col] == (matrix[row-1][col] + gap_score):

            (alignment_left, alignment_top, alignment_position)= self.align(\
                matrix, row-1, col, alignment_left, alignment_top,\
                alignment_position, gap_score, match, mismatch\
            )
        
            alignment_position+= 1  
            alignment_top[alignment_position]= self.buffer[1][row-1].code
            alignment_left[alignment_position]= '-'

        return (alignment_left, alignment_top, alignment_position)

if __name__ == '__main__':

    alignment= Alignment(streams= 2)
    alignment.feed('CGTGAATTCAT', 0)
    alignment.feed('GACTTAC', 1)

    print 'sequences...'
    print [n.code for n in alignment.buffer[0]]
    print [n.code for n in alignment.buffer[1]]
    print


    matrix= alignment.compute_simularity("needleman_wunsch")
    print 'needleman-wunsch matrix...'
    pprint(matrix)
    print

    cols= len(alignment.buffer[0]) # 1 
    rows= len(alignment.buffer[1]) # 1

    alignment_left= []
    alignment_top= []
    for col in range(0, max(cols, rows)):
        alignment_left.append('-')
        alignment_top.append('-')


    #alignment_left= ['-'] * cols 
    #alignment_top= ['-'] * cols
   
    alignment.align(matrix, rows, cols, alignment_left, alignment_top)

    print 'alignments...'
    print alignment_left
    print alignment_top
    print
