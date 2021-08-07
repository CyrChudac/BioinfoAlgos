import os
import sys
from Bio import SeqIO
class DifferentLengthError(Exception):
    """raised when two sequences passed to hamming distance have different length"""
    pass
class HammingDistance:
    def hd_help():
        print("hamming distance help:")
        t = "\t"
        print(t + "tool -hd [sequence1] [s1_number] [sequence2] [s2_number]")
        print(2*t + "both sequence paramaters have to be a path to a fasta file")
        print(2*t + "sequence_numbers are the numbers indicating which sequence to read in given file")
        print()
        print(2*t + "if the sequences length differ, the DifferentLengthError is raised")
    def run(self, ab):
        if len(sys.argv) <= self.ab+2:
            hd_help()