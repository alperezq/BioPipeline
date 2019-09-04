import os
import subprocess


depProc = subprocess.Popen(["BayesTraitsV3", "projects/newTest/BAYESfiles/Ktree.nexus", "projects/newTest/BAYESfiles/OG0000000_mat.txt"], stdin=subprocess.PIPE, text=True)
depCom=("3\n2\nPriorAll exp 10\nStones 100 500\nCor 1\nRun\n")
depProc.stdin.write(depCom)
