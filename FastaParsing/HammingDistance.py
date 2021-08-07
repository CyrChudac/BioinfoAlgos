import os
import sys
import Task
import FastaParsing

class HammingDistance(Task.Task):
    def __init__(self):
        self.name="hamming distance"
        return super(HammingDistance, self).__init__()
    class DifferentLengthError(Exception):
        """raised when two sequences passed to hamming distance have different length"""
        pass
    def hd_help():
        print("hamming distance help:")
        t = "\t"
        print(t + "tool -hd [sequence1] [s1_number] [sequence2] [s2_number]")
        print(2*t + "both sequence paramaters have to be a path to a fasta file")
        print(2*t + "sequence_numbers are the numbers indicating which sequence to read in given file")
        print()
        print(2*t + "if the sequences length differ, the DifferentLengthError is raised")
    def run(self, ab):
        if len(sys.argv) < ab+4:
            hd_help()
        else:
            fp = FastaParsing.FastaParsing()
            s1 = fp.seqInfo(sys.argv[ab+1],sys.argv[ab+2]).seq
            s2 = fp.seqInfo(sys.argv[ab+3],sys.argv[ab+4]).seq
            if len(s1) != len(s2):
                raise DifferentLengthError()
            res = 0
            for i in range(len(s1)):
                if s1[i] != s2[i]:
                    res += 1
            print(res)
        