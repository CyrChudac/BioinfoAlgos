import Task
from Bio import SeqIO

class FastaParsing(Task.Task):
    def __init__(self):
        self.name = "fasta parsing"
    def help(self, pre: str, mult: int):
        print(pre*mult + "[sequence] [sequence_number] ")
        print(pre*(mult+1) + "'sequence' has to be a path to a fasta file")
        print(pre*(mult+1) + "'sequence_number' is the number indicating which sequence to read")
        print()
        print(pre*(mult+1) + "-s --show [d|s|sd|ds]")
        print(pre*(mult+2) + "show sequence and/or it's description in given order")
        print(pre*(mult+1) + "-l --length")
        print(pre*(mult+2) + "displays full sequence length")
        print(pre*(mult+1) + "-ss --subsequence [from] [to]")
        print(pre*(mult+2) + "subsequence from the initial sequence starting at 'from' and ending at 'to'")
        print(pre*(mult+2) + "'from' required, 'to' is not")
    def seqInfo(self, path: str, num: str):
        iter = SeqIO.parse(path,"fasta")
        for i in range(int(num)):
            res = next(iter)
        return res
    def run(self, params):
        if len(params) <= 2:
            self.help()
        else:
            try:
                rec = self.seqInfo(params[1], params[2])
                if params[3] == "-s" or params[3] == "--show":
                    if len(params) == 4:
                        print(rec.description)
                        print(rec.seq)
                    if params[4] == "d":
                        print(rec.description)
                    else:
                        if params[4] == "s":
                            print(rec.seq)
                        else:
                            if params[4] == "sd":
                                print(rec.seq)
                                print(rec.description)
                            else:
                                if params[4] == "ds":
                                    print(rec.description)
                                    print(rec.seq)
                                else:
                                    print("unsupported argument indicating what to show: " + "\"" + params[4] + "\"")
                elif params[3] == "-l" or params[3] == "--length":
                    print(len(rec.seq))
                elif params[3] == "-ss" or params[3] == "-subsequence":
                    if len(params) < 5:
                        self.help()
                    else:
                        if len(params) == 5:
                            print(rec.seq[int(params[4]):])
                        else:
                            print(rec.seq[int(params[4]): int(params[5])])
                else:
                    print("argument " + params[3] + " is not a valid argument:")
                    Task.Task.help(self)
            except StopIteration:
                print("there is no sequence of number " + params[2])
            except ValueError:
                print("argument \"" + params[2] + "\" is not a number")
    
