import os
import subprocess

ProkkaFiles = "projects/halfTest/PROKKAfiles"

for file in ProkkaFiles:
    if os.path.isfile(ProkkaFiles + file):
        print(file)
    else:
        print("False")
