#!/usr/bin/env python3
import pandas as pd
import shutil

#Checks CSV files to be used for bound matrix for missing columns based on FASTAlist
def csvFix(Rfiles, FASTAlist):
    matCSV = pd.read_csv(Rfiles + "matTRF.csv", index_col=0)
    orthoCSV = pd.read_csv(Rfiles + "OrthoGCMatrix.csv", index_col=0)
    groupCSV = pd.read_csv(Rfiles + "distal_GroupMatrix.csv", index_col=0)
    FASTAprefix = []
    for name in FASTAlist:
        temp = name.split(".")
        FASTAprefix.append(temp[0])
    matCols = list(matCSV.columns)
    orthoCols = list(orthoCSV.columns)
    groupCols = list(groupCSV.columns)
    matCompare = (list(set(FASTAprefix).difference(matCols)))
    orthoCompare = (list(set(FASTAprefix).difference(orthoCols)))
    groupCompare = (list(set(FASTAprefix).difference(groupCols)))
    for col in matCompare:
        print("Missing " + col + " in matTRF.csv. Adding, setting values to 0")
        matCSV[col] = int(0)
    for col in orthoCompare:
        print("Missing " + col + " in OrthoGCMatrix.csv. Adding, setting values to 0")
        orthoCSV[col] = int(0)
    for col in groupCompare:
        print("Missing " + col + " in distal_GroupMatrix.csv. Adding, setting values to 0")
        groupCSV[col] = int(0)
    matCSV.to_csv(Rfiles + "matTRF.csv")
    orthoCSV.to_csv(Rfiles + "OrthoGCMatrix.csv")
    groupCSV.to_csv(Rfiles + "distal_GroupMatrix.csv")

#Function for initializing scoary and second round of R modifications
def scoary(pipePath, providedCSV):
    shutil.copyfile(providedCSV, pipePath + "Rfiles/" + providedCSV)
    scorFile = (pipePath + "Rfiles/" + providedCSV)
    if scorFile.lower().endswith(".csv"):
        boundCSV = pd.read_csv(pipePath + "Rfiles/boundMatrix.csv", index_col=0)
        traitCSV = pd.read_csv(scorFile, index_col=0)
        boundCols = list(boundCSV.columns)
        traitRows = list(traitCSV.index)
        boundTraitCompare = list(set(boundCols) - set(traitRows))
        if 'V1' in boundTraitCompare:
            boundTraitCompare.remove('V1')
        if len(boundTraitCompare) is 0:
            return scorFile;
        else:
            print("Columns of boundMatrix do not match rows of provided CSV, unable to run Scoary")
    else:
        print("Provided csv does not end with csv extension, unable to proceed.")
