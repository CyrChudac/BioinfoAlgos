import Task
import MsaParsing
from Bio import AlignIO
import math

class PositionConservation(Task.Task):
    threshold = 0.17
    window = 5
    def __init__(self):
        self.name = "msa position conservation"
    def help(self, pre: str, mult: int):
        print(pre*mult + "[alignmnet] [scoring_matrix]")
        print(pre*(mult+1) + "'alignmnet' is a path to the alignment file")
        print(pre*(mult+1) + "'scoring_matrix' is a path to the scoring matrix file")
        print()
        print(pre*(mult+1) + "-f --first [n]")
        print(pre*(mult+2) + "first n best scoring columns will have their number and score displayed")
        print(pre*(mult+1) + "-ic --is-conservated [n]")
        print(pre*(mult+2) + "returns true, if the average conservation rate is below threshold")
        print(pre*(mult+2) + "among the n-th residue")
        print()
        print(pre*(mult+2) + "-w --window [x]")
        print(pre*(mult+3) + "x specifies how many residues to each side from the n-th residue")
        print(pre*(mult+3) + "will be considered; default value is " + str(PositionConservation.window))
        print(pre*(mult+2) + "-t --threshold [x]")
        print(pre*(mult+3) + "x specifies the threshold under which the average has to be to pass")
        print(pre*(mult+3) + "default value is " + str(PositionConservation.threshold))
    def run(self, params):
        if len(params) < 3:
            print("wrong number of parameters:")
            Task.Task.help(self)
        else:
            if not Task.isFile(params[0]) or not Task.isFile(params[1]):
                print("one of given paths is not a path")
                return None
            parser = MsaParsing.MSAParser()
            ali = parser.run([params[0],"-s"]).val
            result = []
            for i in range(0, len(ali[0].seq)):
                sum = parser.run([params[0],"-sop",params[1],str(i)]).val
                result.append((i, sum))
            result.sort(key=lambda t: t[1])
            #removing the negative numbers, if they exist
            if result[0][1] < 0:
                lowest = result[0][1]
                for i in range(0, len(result)):
                    result[i] = (result[i][0], result[i][1] - lowest)

            #n*(n+1)/2, but n is actually (n-1) here, so (n-1)*(n-1+1)/2
            additions = len(ali)*(len(ali)-1)/2
            for i in range(0, len(result)):
                result[i] = (result[i][0], result[i][1]/additions)
            if params[2][0:2] == "-f" or params[2] == "--first":
                if params[2][0:2] == "--" or len(params[2]) == 2:
                    if len(params) == 4:
                        n = int(params[3])
                    else:
                        print("wrong number of parameters:")
                        Task.Task.help(self)
                        return None
                else:
                    n = int(params[2][2:])
                return Task.Result(result[0:n], Task.Result.printList_advanced(self.show))
            elif params[2][0:3] == "-ic" or params[2] == "--is-conservated":
                if len(params) < 4:
                    print("missing number indicating, which position to test")
                    return None
                n = int(params[3])
                result.sort(key=lambda t: t[0])
                t = PositionConservation.threshold
                w = PositionConservation.window
                i = 4
                while i < len(params):
                    if params[i] == "-w" or params[i] == "--window":
                        w = int(params[i+1])
                        i += 1
                    elif params[i] == "-t" or params[i] == "--threshold":
                        t = float(params[i+1])
                        i += 1
                    i += 1
                start = max(0, n - w)
                end = min(len(result), n + w + 1)
                sum = 0
                for i in range(start, end):
                    sum += result[i][1]
                sum /= (end-start)
                return Task.Result((sum,t), PositionConservation.showIc)
            else:
                print("argument " + params[2] + " is not a valid argument:")
                Task.Task.help(self)
    def show(self,x):
        return "column " + str(x[0]) + ":\t" + str(x[1])
    def showIc(x):
        sum, t = x
        print("(" + str(sum) + " < " + str(t) + ") = " + str(sum < t))