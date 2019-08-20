#!/usr/bin/env python3
import pandas as pd
import multiprocessing as mp
import subprocess
from itertools import repeat
import shutil
import os
import re

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
    shutil.copyfile(providedCSV, pipePath + "Rfiles/providedCSV.csv")
    scorFile = (pipePath + "Rfiles/providedCSV.csv")
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

#Parsing of BayesGenerator.csv, then shuttling of data to Bayes
def bayesPool(pipePath):
    bayesData = pd.read_csv(pipePath + "BAYESfiles/BayesGenerator.csv", index_col=0)
    with open(pipePath + "Results/bayesResults.txt", 'w') as bayesOut:
        bayesOut.writelines("Gene\tDependent\tIndependent\n")
    bayesList = (bayesData.columns)
    pool = mp.Pool(processes=20,)
    pool.starmap(bayesCall, zip(bayesList, repeat(pipePath)))
    pool.close()
    pool.join()

#Helper function for pool process of BAYESfiles
def bayesCall(iterableColumns, path):
    nexus = (path + "BAYESfiles/Ktree.nexus")
    bayesData = pd.read_csv(path + "BAYESfiles/BayesGenerator.csv", index_col=0)
    traitData = pd.read_csv(path + "Rfiles/providedCSV.csv", index_col=0)
    df = bayesData[[iterableColumns]].join(traitData)
    dfName = (path + "BAYESfiles/" + iterableColumns + "_mat.txt")
    df.to_csv(dfName, header=False, sep= "\t")
    os.system("BayesTraitsV3 " + nexus + " " + dfName + " < addScripts/addResources/Comm_dependant")
    os.rename(dfName + ".Stones.txt", path + "BAYESfiles/" + iterableColumns + ".dependant")
    os.system("BayesTraitsV3 " + nexus + " " + dfName + " < addScripts/addResources/Comm_independant")
    os.rename(dfName+ ".Stones.txt", path + "BAYESfiles/" + iterableColumns + ".independant")
    with open(path + "BAYESfiles/" + iterableColumns + ".dependant", 'r') as depFile:
        content = depFile.readlines()
        for line in content:
            if "Log marginal likelihood" in line:
                temp = line.split()
                tempDep = temp[3]
    with open(path + "BAYESfiles/" + iterableColumns + ".independant", 'r') as indepFile:
        content = indepFile.readlines()
        for line in content:
            if "Log marginal likelihood" in line:
                temp = line.split()
                tempInd = temp[3]
    try:
        tempDep
    except NameError:
        print("Missing value for " + iterableColumns + " in dependant file.")
        dependantVal = ' '
    else:
        dependantVal = tempDep
    try:
        tempInd
    except NameError:
        print("Missing value for " + iterableColumns + " in independant file.")
        independantVal = ' '
    else:
        independantVal = tempInd
    with open(path + "Results/bayesResults.txt", 'a') as bayesOut:
        bayesOut.writelines(iterableColumns + "\t" + dependantVal + "\t" + independantVal + "\n")
    for item in os.listdir(path + "BAYESfiles/"):
        if iterableColumns in item:
            os.remove(path + "BAYESfiles/" + item)
