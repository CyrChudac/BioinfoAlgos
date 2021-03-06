import Task
import FastaParsing

class HammingDistance(Task.Task):
    def __init__(self):
        self.name="hamming distance"
        return super(HammingDistance, self).__init__()
    def help(self, pre:str, mult: int):
        print(pre*mult + "[sequence1] [s1_number] [sequence2] [s2_number]")
        print(pre*(mult+1) + "both sequence paramaters have to be a path to a fasta file")
        print(pre*(mult+1) + "sequence_numbers are the numbers indicating which sequence to read in given file")
    def run(self, params):
        if len(params) < 4:
            print("wrong number of parameters:")
            Task.Task.help(self)
        else:
            if not Task.isFile(params[0]) or not Task.isFile(params[2]):
                print("one of given paths is not a path")
                return None
            fp = FastaParsing.FastaParsing()
            s1 = fp.seqInfo(params[0],params[1]).seq
            s2 = fp.seqInfo(params[2],params[3]).seq
            if len(s1) != len(s2):
                print("hammming distance can only be computed for sequences of the same length(" 
                      + str(len(s1)) + " != " + str(len(s2)) + ")")
                return None
            res = 0
            for i in range(len(s1)):
                if s1[i] != s2[i]:
                    res += 1
            return Task.Result(res)
        