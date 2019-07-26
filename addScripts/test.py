#!/usr/bin/python
import os
import pandas as pd

def concatNuc(rvdFiles, results):
    nucList = [file for file in os.listdir(rvdFiles) if file.endswith.csv]
    df = pd.concat(map(pd.read_csv, nucList))
    print(df)

rvdFiles = "projects/idunsTEST/RVDfiles/"
results = "projects/idunsTEST/Results/"

concatNuc(rvdFiles, results)
