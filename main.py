import sys
import FastaParsing
import HammingDistance
import EditDistance
import PdbParsing
import MsaParsing
import positionConservation
import MsaStructure
import help
import Task

dic = {
    "-fp": FastaParsing.FastaParsing(),
    "-hd": HammingDistance.HammingDistance(),
    "-ed": EditDistance.EditDistance(),
    "-pdb": PdbParsing.PdbParser(),
    "-msap": MsaParsing.MSAParser(),
    "-pc": positionConservation.PositionConservation(),
    "-msas": MsaStructure.MSAStructure()
}
longer = {
    "--fasta-parsing": "-fp",
    "--hamming-distance": "-hd",
    "--edit-distance": "-ed",
    "--pdb-parsing": "-pdb",
    "--multiseq-align-parsing": "-msap",
    "--position-conservation": "-pc",
    "--multiseq-align-structure": "-msas"}
try:

    if sys.argv[1][0:2] == "--":
        task = dic[longer[sys.argv[1]]]
    else:
        task = dic[sys.argv[1]]
except (KeyError, IndexError):
    task = help.Help(dic, longer)
    
result = task.run(sys.argv[2:])

"""
try:
    result = task.run(sys.argv[2:])
except Exception as err:
    print("!!!!!!!!!Erorr!!!!!!!!!")
    print(sys.argv)
    print(err)
    print("press any key to continue...")
    sys.stdin.readline()
"""
if result != None:
    result.display()
print("press any key to continue...")
sys.stdin.readline()