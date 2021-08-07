import os
import sys
import Task
import FastaParsing

class EditDistance(Task.Task):
    def __init__(self):
        self.name="edit distance"
        self.curr_results = []
        return super(EditDistance, self).__init__()
    def hd_help():
        print("edit distance help:")
        t = "\t"
        print(t + "tool -ed [sequence1] [s1_number] [sequence2] [s2_number] [-a|-m|-ma]")
        print(2*t + "both sequence paramaters have to be a path to a fasta file")
        print(2*t + "sequence_numbers are the numbers indicating which sequence to read in given file")
        print(2*t + "if optional parameter -a is used, all alignments with the best score with be displayed too")
        print(2*t + "if optional parameter -m is used, the dynamic matrix will be displayed")
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

    def run(self, ab):
        if len(sys.argv) < ab+4:
            hd_help()
        else:
            fp = FastaParsing.FastaParsing()
            s1 = fp.seqInfo(sys.argv[ab+1],sys.argv[ab+2]).seq
            s2 = fp.seqInfo(sys.argv[ab+3],sys.argv[ab+4]).seq
            m = self.getMatrix(s1, s2)
            ed = m[len(s1)-1][len(s2)-1]
            print("edit distance is:")
            print(ed)
            if len(sys.argv) == ab+6:
                if sys.argv[ab+5] == "-m" or sys.argv[ab+5] == "-ma":
                    print()
                    print("with matrix:")
                    for i in m:
                        print(i)
                if sys.argv[ab+5] == "-a" or sys.argv[ab+5] == "-ma":
                    cr = self.backtrack(s1,s2,m)
                    print()
                    print("with alignments:")
                    for x in cr:
                        print(x[0])
                        print(x[1])
                        print()


        