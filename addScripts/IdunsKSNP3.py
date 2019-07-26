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
def ksnpParse(scoaryDir, rDir, providedCSV, disDir, trfDir, orthoDir, kDir, resultsDIR, rvdDir):
    for file in os.listdir(scoaryDIR):
        if file.endswith(".csv"):
            scoaryCSV = scoaryDir +file
    repeatsCSV = rDir + "RepeatNames.csv"
    boundFile = rDir + "boundMatrix.csv"
    treeFile = kDir + "kSNP3_results/tree.ML.tre"
    parsedTRF = trfDir + "trfParsed.txt"
    distalGroups = disDir + "Outputs/disTALOut.TALgroups.csv"
    comboFile = disDir + "rvdCombo.FASTA"
    concatenated = concatNuc(rvdDir, rDir)
    subprocess.Popen(["Rscript", "addScripts/IdunsScoaryR.R", scoaryCSV, repeatsCSV, boundFile, treeFile, providedCSV, resultsDIR, parsedTRF, distalGroups, comboFile, concatenated, faaFile], close_fds=True).communicate()[0]


def concatNuc(rvdFiles, results):
    nucList = [file for file in os.listdir(rvdFiles) if file.endswith.csv]
    df = pd.concat(map(pd.read_csv, nucList))
