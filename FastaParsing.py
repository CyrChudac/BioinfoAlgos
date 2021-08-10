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
        print(pre*(mult+1) + "-s --show [-d|-s|-sd|-ds]")
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
    def descSeq(rec):
        print(rec.description)
        print(rec.seq)
    def seq(rec):
        print(rec.seq)
    def desc(rec):
        print(rec.description)
    def seqDesc(rec):
        print(rec.seq)
        print(rec.description)
    def run(self, params):
        if len(params) < 3:
            print("wrong number of parameters:")
            Task.Task.help(self)
        else:
            try:
                if not Task.isFile(params[0]):
                    print("one of given paths is not a path")
                    return None
                rec = self.seqInfo(params[0], params[1])
                if params[2] == "-s" or params[2] == "--show":
                    if len(params) == 3:
                        deleg = FastaParsing.descSeq
                    elif params[3] == "-d":
                        deleg = FastaParsing.desc
                    elif params[3] == "-s":
                            deleg = FastaParsing.seq
                    elif params[3] == "-sd":
                            deleg = FastaParsing.seqDesc
                    elif params[3] == "-ds":
                            deleg = FastaParsing.descSeq
                    else:
                        print("unsupported argument indicating what to show: " + "\"" + params[3] + "\"")
                        return None
                    return Task.Result(rec, deleg)
                elif params[2] == "-l" or params[2] == "--length":
                    return Task.Result(len(rec.seq))
                elif params[2] == "-ss" or params[2] == "-subsequence":
                    if len(params) < 4:
                        self.help()
                        return None
                    else:
                        if len(params) == 4:
                            return Task.Result(rec.seq[int(params[3]):])
                        else:
                            return Task.Result(rec.seq[int(params[3]): int(params[4])])
                else:
                    print("argument " + params[2] + " is not a valid argument:")
                    Task.Task.help(self)
            except StopIteration:
                print("there is no sequence of number " + params[1])
            except ValueError:
                print("argument \"" + params[1] + "\" is not a number")
    
