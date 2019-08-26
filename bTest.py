import os
import subprocess


P = subprocess.Popen(["BayesTraitsV3", "projects/halfTest/BAYESfiles/Ktree.nexus", \
"projects/halfTest/BAYESfiles/OG0004270_mat.txt"], \
stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
#P.stdin.write("3")
#P.stdin.write("2")
P.stdin.write("3\n2\nPriorAll exp 10\nStones 100 500\nCor 1\nRun\n")
#P.stdin.write("Stones 100 500\n")
#P.stdin.write("Cor 1\n")
#P.stdin.write("Run\n")
