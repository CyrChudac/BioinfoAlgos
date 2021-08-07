import sys
import FastaParsing
import HammingDistance
import EditDistance
import Task

task = Task.Task()
if sys.argv[1] == "-fp":
    task = FastaParsing.FastaParsing()
if sys.argv[1] == "-hd":
    task = HammingDistance.HammingDistance()
if sys.argv[1] == "-ed":
    task = EditDistance.EditDistance()
task.run(1)