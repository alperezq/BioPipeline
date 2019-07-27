#!/usr/bin/env python3
import os
import shutil

def concatFaa(prokkaDir, results):
    faaList = [file for file in os.listdir(prokkaDir) if file.lower().endswith(".faa")]
    faaList.sort()
    for file in faaList:
        prefix = file.split(".")
        with open(prokkaDir + file, 'r') as item:
            data = item.read()
        data = data.replace('>', '>' + prefix[0] + ' ')
        with open(prokkaDir + file, 'w') as item:
            item.write(data)
    with open(results + "faaConcatenated.faa", 'wb') as outFile:
        for file in faaList:
            with open(prokkaDir + file, 'rb') as inFile:
                shutil.copyfileobj(inFile, outFile)


prokkaDir = "projects/idunsTEST/PROKKAfiles/FAAsTEST/"
results = "projects/idunsTEST/Results/"
concatFaa(prokkaDir, results)
