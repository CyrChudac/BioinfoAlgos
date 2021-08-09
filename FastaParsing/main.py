import sys
import FastaParsing
import HammingDistance
import EditDistance
import PdbParser
import help
import Task

dic = {
    "-fp": FastaParsing.FastaParsing(),
    "-hd": HammingDistance.HammingDistance(),
    "-ed": EditDistance.EditDistance(),
    "-pdb": PdbParser.PDBParsing(),
    }
dic["-h"] = help.Help(dic)
try:
    task = dic[sys.argv[1]]
except KeyError:
    task = dic["-h"]

task.run(1)