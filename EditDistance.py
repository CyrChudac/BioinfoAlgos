import sys
import Task
import FastaParsing

class EditDistance(Task.Task):
    def __init__(self):
        self.name="edit distance"
        self.curr_results = []
        return super(EditDistance, self).__init__()
    def help(self, pre: str, mult: int):
        print(pre*mult + "[sequence1] [s1_number] [sequence2] [s2_number] ")
        print(pre*(mult+1) + "both sequence paramaters have to be a path to a fasta file")
        print(pre*(mult+1) + "sequence_numbers are the numbers indicating which sequence to read in given file")
        print()
        print(pre*(mult+1) + "-a")
        print(pre*(mult+2) + "display alignments")
        print(pre*(mult+1) + "-m")
        print(pre*(mult+2) + "display the dynamic matrix used for the calculation")
        print(pre*(mult+1) + "-ma")
        print(pre*(mult+2) + "is used if both alignments and matrix should be shown")
    def getMatrix(self, s1, s2):
        h = len(s1) + 1
        w = len(s2) + 1
        m = [[0] * w for i in range(h)]
        for i in range(1,h):
            m[i][0] = i
        for j in range(1,w):
            m[0][j] = j
        for i in range(1,h):
            for j in range(1,w):
                r = sys.maxsize
                r = min(r, m[i][j-1] + 1)
                r = min(r, m[i-1][j] + 1)
                if s1[i - 1] == s2[j - 1]:
                    r = min(r, m[i-1][j-1])
                else:
                    r = min(r, m[i-1][j-1] + 1)
                m[i][j] = r
        return m
    def backtrackRecur(self, s1, i1, r1, s2, i2, r2, m):
        if i1 == 0 and i2 == 0:
            self.curr_results.append([r1,r2])
        else:
            if i1 > 0:
                if m[i1][i2] == m[i1-1][i2] + 1:
                    self.backtrackRecur(s1,i1-1,s1[i1-1] + r1,s2,i2,"-" + r2,m)
            if i2 > 0:
                if m[i1][i2] == m[i1][i2-1] + 1:
                    self.backtrackRecur(s1,i1,"-" + r1,s2,i2-1,s2[i2-1] + r2,m)
            if i1 > 0 and i2 > 0:
                if s1[i1-1] == s2[i2-1]:
                    if m[i1][i2] == m[i1-1][i2-1]:
                        self.backtrackRecur(s1,i1-1,s1[i1-1] + r1,s2,i2-1,s2[i2-1] + r2,m)
                else:
                    if m[i1][i2] == m[i1-1][i2-1] + 1:
                        self.backtrackRecur(s1,i1-1,s1[i1-1] + r1,s2,i2-1,s2[i2-1] + r2,m)
    def backtrack(self, s1, s2, m):
        self.curr_results.clear()
        r1= ""
        r2 = ""
        self.backtrackRecur(s1,len(s1),"",s2,len(s2),"",m)
        return self.curr_results
    
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
            m = self.getMatrix(s1, s2)
            ed = m[len(s1)-1][len(s2)-1]
            if len(params) == 5:
                if params[4] == "-m":
                    return Task.Result((ed, m, None), EditDistance.m)
                elif params[4] == "-a":
                    return Task.Result((ed, m, None), EditDistance.a)
                elif params[4] == "-ma":
                    return Task.Result((ed, m, None), EditDistance.ma)
                else:
                    print("argument " + params[4] + " is not a valid argument:")
                    Task.Task.help(self)
            else:
                return Task.Result((ed, m, None), EditDistance.alignShow)

    def distanceShow(trip):
        ed, _, _ = trip
        print("edit distance is:")
        print(ed)
    def matrixShow(trip):
        ed, m, cr = trip
        print()
        print("with matrix:")
        for i in m:
            print(i)
    def alignShow(trip):
        ed, m, cr = trip
        print()
        print("with alignments:")
        for x in cr:
            print(x[0])
            print(x[1])
            print()
    def m(trip):
        EditDistance.distanceShow(trip)
        EditDistance.matrixShow(trip)
    def ma(trip):
        EditDistance.distanceShow(trip)
        EditDistance.matrixShow(trip)
        EditDistance.alignShow(trip)
    def a(trip):
        EditDistance.distanceShow(trip)
        EditDistance.alignShow(trip)