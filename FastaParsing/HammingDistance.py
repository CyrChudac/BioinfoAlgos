import Task
import FastaParsing

class HammingDistance(Task.Task):
    def __init__(self):
        self.name="hamming distance"
        return super(HammingDistance, self).__init__()
    class DifferentLengthError(Exception):
        """raised when two sequences passed to hamming distance have different length"""
        pass
    def help(self, pre:str, mult: int):
        print(pre*mult + "[sequence1] [s1_number] [sequence2] [s2_number]")
        print(pre*(mult+1) + "both sequence paramaters have to be a path to a fasta file")
        print(pre*(mult+1) + "sequence_numbers are the numbers indicating which sequence to read in given file")
        print()
        print(pre*(mult+1) + "if the sequences length differ, the DifferentLengthError is raised")
    def run(self, params):
        if len(params) < 4:
            self.help()
        else:
            fp = FastaParsing.FastaParsing()
            s1 = fp.seqInfo(params[1],params[2]).seq
            s2 = fp.seqInfo(params[3],params[4]).seq
            if len(s1) != len(s2):
                raise DifferentLengthError()
            res = 0
            for i in range(len(s1)):
                if s1[i] != s2[i]:
                    res += 1
            print(res)
        