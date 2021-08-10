import Task
from Bio import AlignIO

class MSAParser(Task.Task):
    class ScoringMatrix:
        def __init__(self, defaultOk = 0, defaultNot = 1):
            self.matrix = {"default": {"ok": defaultOk, "not": defaultNot}}
        
        def fromFile(path: str):
            f = open(path)
            s = f.readline()
            if "\t" not in s:
                matrix = MSAParser.ScoringMatrix(int(s.split("|")[0]), int(s.split("|")[1]))
                s = f.readline()
            else:
                matrix = ScoringMatrix()
            s = s.split("\t")
            for x in f:
                data = x.split("\t")
                data = data[0:len(data)-1]
                for i in range(1,len(data)):
                    matrix[data[0],s[i]] = float(data[i])
            return matrix
        def __setitem__(self, tup, val):
            key1, key2 = tup
            if not key1 in self.matrix.keys():
                self.matrix[key1]={}
            self.matrix[key1][key2] = val
            if not key2 in self.matrix.keys():
                self.matrix[key2]={}
            self.matrix[key2][key1] = val
        def __getitem__(self, tup):
            key1, key2 = tup
            if key1 in self.matrix.keys() and key2 in self.matrix[key1].keys():
                return self.matrix[key1][key2]
            else:
                if key1 == key2:
                    return self.matrix["default"]["ok"]
                else:
                    return self.matrix["default"]["not"]
    def __init__(self):
        self.name = "msa parsing"
    def help(self, pre: str, mult: int):
        print(pre*mult + "[file]")
        print(pre*(mult+1) + "'file' is a path to the alignment file")
        print()
        print(pre*(mult+1) + "-s --show")
        print(pre*(mult+2) + "displays the whole alignment")
        print(pre*(mult+1) + "-sq --sequence [-id|-n] [indicator]")
        print(pre*(mult+2) + "displays sequence given by id or N-th sequence in the alignment")
        print(pre*(mult+2) + "indicator is the value of id respectively of N")
        print(pre*(mult+1) + "-c --column [n]")
        print(pre*(mult+2) + "displays a whole column of the alignment")
        print(pre*(mult+2) + "'n' is a zero-based index of the column")
        print(pre*(mult+1) + "-sop --sum-of-pairs [scoring_matrix] [n]")
        print(pre*(mult+2) + "displays the value of sum of pairs")
        print(pre*(mult+2) + "'scoring matrix' - optional - path to a file with a scoring matrix")
        print(pre*(mult+2) + "'n' - optional - if given, only sum of pairs of N-th column will be computed")
    def findSeq(self, ali, id):
        for s in ali:
            if s.id == id:
                return s
    def column(self, ali, c):
        result = []
        for s in ali:
            result.append(s.seq[c])
        return result
    def run(self, params):
        if len(params) < 2:
            Task.task.help(self)
        else:
            align = AlignIO.read(params[0],"clustal")
            if params[1] == "-s" or params[1] == "--show":
                return Task.Result(align)
            elif params[1] == "-sq" or params[1] == "--sequence":
                if params[2][0:2] == "-n":
                    if len(params[2]) == 2:
                        if len(params) > 3:
                            result = align[int(params[3])].seq
                        else:
                            print("sequence number not given")
                            return None
                    else:
                        result = align[int(params[2][2:])].seq
                if params[2] == "-id":
                    if len(params) > 3:
                        result = self.findSeq(align, params[3]).seq
                    else:
                        print("sequence number not given")
                return Task.Result(result)
            elif params[1][0:2] == "-c" or params[1] == "--column":
                if len(params) == 3:
                    n = int(params[2])
                elif params[1][0:2] == "-c" and len(params[1]) != 2:
                    n = int(params[1][2:])
                else:
                    print("wrong number of parameters")
                    return None
                arr = self.column(align, n)
                for a in arr:
                    print(a)
                return Task.Result(arr, Task.Result.printList)
            elif params[1] == "-sop" or params[1] == "--sum-of-pairs":
                start=0
                end=len(align[0])
                if len(params) == 2:
                    matrix = self.ScoringMatrix()
                else:
                    matrix = MSAParser.ScoringMatrix.fromFile(params[2])
                    if len(params) > 3:
                        start = int(params[3])
                        end = start + 1
                sum = 0
                for c in range(start, end):
                    col = self.column(align, c)
                    for i in range(0, len(col)):
                        for j in range(0, i):
                            sum += matrix[col[i], col[j]]
                return Task.Result(sum)
            else:
                print("argument " + params[1] + " is not a valid argument:")
                Task.Task.help(self)
