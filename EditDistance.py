import Task
import FastaParsing

class EditDistance(Task.Task):
    def __init__(self):
        self.name="edit distance"
        return super(EditDistance, self).__init__()
    def help(self, pre: str, mult: int):
        print(pre*mult + "[sequence1] [s1_number] [sequence2] [s2_number] ")
        print(pre*(mult+1) + "both sequence paramaters have to be a path to a fasta file")
        print(pre*(mult+1) + "sequence_numbers are the numbers indicating which sequence to read in given file")
        print()
        print(pre*(mult+1) + "-m")
        print(pre*(mult+2) + "display the dynamic matrix used for the calculation; no backtrack computed")
        print(pre*(mult+1) + "-a")
        print(pre*(mult+2) + "display alignments")
        print(pre*(mult+1) + "-ma")
        print(pre*(mult+2) + "is used if both alignments and matrix should be shown")
        print(pre*(mult+1) + "--only [x]")
        print(pre*(mult+2) + "if specified with -a or -ma, only x backtrack results will be computed")
    def getMatrix(self, s1, s2):
        h = len(s1) + 1
        w = len(s2) + 1
        m = [[[0,[[0,0]]]] * w for i in range(h)]
        for i in range(1,h):
            m[i][0] = [i,[[-1,0]]]
        for j in range(1,w):
            m[0][j] = [j,[[0,-1]]]
        for i in range(1,h):
            for j in range(1,w):
                r = [m[i][j-1][0] + 1,[[0,-1]]]
                if  m[i-1][j][0] + 1 < r[0]:
                    r = [m[i-1][j][0] + 1,[[-1,0]]]
                elif m[i-1][j][0] + 1 == r[0]:
                    r[1].append([-1,0])

                if s1[i - 1] == s2[j - 1]:
                    k = m[i-1][j-1][0]
                else:
                    k = m[i-1][j-1][0] + 1
                
                if k < r[0]:
                    r = [k,[[-1,-1]]]
                elif k == r[0]:
                    r[1].append([-1,-1])
                m[i][j] = r
        return m
    def fromPath(self, path, s1, s2):
        r1 = ""
        r2 = ""
        l1 = len(s1)
        l2 = len(s2)
        for i in path[1:]:
            if i[0] == l1:
                r1 = "-" + r1
            else:
                r1 = s1[i[0]] + r1
            if i[1] == l2:
                r2 = "-" + r2
            else:
                r2 = s2[i[1]] + r2
            l1 = i[0]
            l2 = i[1]
        return [r1,r2]
    def backtrack(self, s1, s2, m, only):
        curr_results = []
        print("# matrix computed")
        stack = [(len(s1),len(s2),0)]
        while(len(stack) > 0):
            i1, i2, i3 = stack.pop()
            if i1 == 0 and i2 == 0:
                stack.append([i1,i2])
                curr_results.append(self.fromPath(stack, s1, s2))
                stack.pop()
                if len(curr_results) == only:
                    return curr_results
                if only >= 2000 and len(curr_results) % int(only/10) == 0:
                    print("# " + str(int(len(curr_results) * 100 / only) + 1) + "%")
            elif i3 < len(m[i1][i2][1]):
                stack.append((i1,i2,i3+1))
                coo = m[i1][i2][1][i3]
                stack.append((i1+coo[0],i2+coo[1],0))
        return curr_results
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
            numsM = [[j[0] for j in i] for i in m]
            ed = m[len(s1)-1][len(s2)-1][0]
            if len(params) >= 5:
                if params[4] == "-m":
                    return Task.Result((ed, numsM, None),
                       lambda x: self.distanceShow(x).matrixShow(x))
                else:
                    if len(params) >= 7 and params[5] == "--only":
                        cr = self.backtrack(s1, s2, m, int(params[6]))
                    else:
                        cr = self.backtrack(s1, s2, m, -1)
                    if params[4] == "-a":
                        return Task.Result((ed, numsM, cr), 
                           lambda x: self.distanceShow(x).alignShow(x))
                    elif params[4] == "-ma":
                        return Task.Result((ed, numsM, cr), 
                           lambda x: self.distanceShow(x).matrixShow(x).alignShow(x))
                    else:
                        print("argument " + params[4] + " is not a valid argument:")
                        Task.Task.help(self)
            else:
                return Task.Result((ed, numsM, None), self.distanceShow)

    def distanceShow(self, trip):
        ed, _, _ = trip
        print("edit distance is:")
        print(ed)
        return self
    def matrixShow(self, trip):
        ed, m, cr = trip
        print()
        print("with matrix:")
        l = len(m[0])
        lens = [0] * l
        for i in m:
            for j in range(l):
                lens[j] = max(lens[j], len(str(i[j])))
        for i in m:
            print("[", end='')
            for j in range(l):
                if j != 0:
                    print(",", end='')
                print(" "*(lens[j]-len(str(i[j]))) + str(i[j]), end='')
            print("]")
        return self
    def alignShow(self, trip):
        ed, m, cr = trip
        print()
        print("with alignments:")
        for x in cr:
            print("> 1")
            Task.Result.printSeq(x[0])
            print("> 2")
            Task.Result.printSeq(x[1])
            print()
        return self