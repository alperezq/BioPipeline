 import subprocess
import sys
import shutil
import os
import re


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


#Function for parsing of KSNP3 Results
def ksnpParse(scoaryDIR, boundFile, treeFile, providedCSV, resultsDIR, parsedTRF, repeatsCSV, distalGroups, comboFile):
    for file in os.listdir(scoaryDIR):
        if file.endswith(".csv"):
            scoaryCSV = scoaryDir +file
    subprocess.Popen(["Rscript", "addScripts/IdunsScoaryR.R", scoaryCSV, boundFile, treeFile, providedCSV, resultsDIR, parsedTRF, repeatsCSV, distalGroups, comboFile, concatenated, faaFile], close_fds=True).communicate()[0]
