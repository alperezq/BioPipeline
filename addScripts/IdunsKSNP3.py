#!/usr/bin/env python3
import subprocess
import sys
import shutil
import os
import re
import pandas as pd


sys.dont_write_bytecode = True


#Function for use of KSNP3, appropriate parsing, and movement of files
def ksnpCall(faPath, ksnpPath, ksnpList, ksnpCpus):
    ksnpGenomes = ksnpPath + "ksnpGenomes.txt"
    with open(ksnpGenomes, "w") as outFile:
        for file in ksnpList:
            prefixKSNP = file.split(".")
            outFile.write(os.path.abspath(faPath + file) + "\t" + prefixKSNP[0] + "\n")
    subprocess.Popen(["MakeFasta", ksnpGenomes, ksnpPath + "ForKchoser"], close_fds=True).communicate()[0]
    subprocess.Popen(["Kchooser", ksnpPath + "ForKchoser"], close_fds=True).communicate()[0]
    shutil.move("Kchooser.report", ksnpPath + "Kchooser.report")
    with open(ksnpPath + "Kchooser.report", "r") as kOut:
        content = kOut.readlines()
        for line in content:
            if re.match(r'The optimum value of K is (.*)\.', line):
                tempInt = re.findall(r'The optimum value of K is (.*)\.', line)
                ksnpInt = tempInt[0]
    subprocess.Popen(["kSNP3", "-in", ksnpGenomes, "-k", str(ksnpInt), "-outdir", ksnpPath + "kSNP3_results", "-ML", "-CPU", str(ksnpCpus)], close_fds=True).communicate()[0]
    shutil.move("fasta_list", ksnpPath + "fasta_list")


#Function for routing needed data to R for parsing of various data sets
def ksnpParse(scoaryDir, rDir, providedCSV, disDir, trfDir, orthoDir, kDir, resultsDir, rvdDir, faaDir):
    for file in os.listdir(scoaryDir):
        if file.endswith(".csv"):
            scoaryCSV = scoaryDir +file
    repeatsCSV = rDir + "RepeatNames.csv"
    boundFile = rDir + "boundMatrix.csv"
    distalGroups = disDir + "Outputs/disTALOut.TALgroups.csv"
    parsedTRF = trfDir + "trfParsed.txt"
    orthogroupsTXT = orthoDir + "Orthogroups.txt"
    treeFile = kDir + "kSNP3_results/tree.ML.tre"
    comboFile = disDir + "rvdCombo.FASTA"
    concatenated = concatNuc(rvdDir, resultsDir)
    faaFile = concatFaa(faaDir, resultsDir)
    subprocess.Popen(["Rscript", "addScripts/IdunsSecondR.R", scoaryCSV, repeatsCSV, boundFile, providedCSV, distalGroups, parsedTRF, orthogroupsTXT, treeFile, comboFile, resultsDir, concatenated, faaFile], close_fds=True).communicate()[0]

#concatenates rvdNuc files to single CSV
def concatNuc(rvdFiles, results):
    nucList = []
    for file in os.listdir(rvdFiles):
        if file.endswith(".csv"):
            nucList.append(rvdFiles + file)
    nucList.sort()
    combined = pd.concat([pd.read_csv(file) for file in nucList])
    combined.to_csv(results + "rvdNucs.csv", index=False)
    return (results + "rvdNucs.csv")

#Creates concatenated FAA file, adding file names to start of appropriate lines
def concatFaa(faaDir, results):
    faaList = [file for file in os.listdir(faaDir) if file.lower().endswith(".faa")]
    faaList.sort()
    for file in faaList:
        prefix = file.split(".")
        with open(faaDir + file, 'r') as item:
            data = item.read()
        data = data.replace('>', '>' + prefix[0] + ' ')
        with open(faaDir + file, 'w') as item:
            item.write(data)
    with open(results + "faaConcatenated.faa", 'wb') as outFile:
        for file in faaList:
            with open(faaDir + file, 'rb') as inFile:
                shutil.copyfileobj(inFile, outFile)
    return(results + "faaConcatenated.faa")
