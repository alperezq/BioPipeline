#!/usr/bin/env python3
import subprocess
import sys
import shutil
import os
import re
import pandas as pd
import multiprocessing as mp
from itertools import repeat

sys.dont_write_bytecode = True

#Parsing of TRF files to single text file
def trfParse(trfSrc, fileList):
    for datFile in os.listdir(os.getcwd()):
        for endFile in fileList:
            if endFile in datFile:
                shutil.move(os.getcwd() + "/" + datFile, trfSrc + datFile)
    with open(trfSrc + "trfParsed.txt", "w") as outFile:
        for file in os.listdir(trfSrc):
            if file.endswith(".dat"):
                prefixTRF = file.split(".")[0]
                with open(trfSrc + file) as inFile:
                    for line in inFile.readlines():
                        if (len(line) >= 45) and ("Sequence" not in line):
                            entries = line.split()
                            entry = (prefixTRF + " " + entries[0] + " " + entries[1] + " " + entries[2] + " " + entries[3] + " " + entries[4] + " " + entries[13] + " " + entries[14] + "\n")
                            outFile.write(entry)


#Run FASTA files through prokka, parse FAA'srun results through Ortho, move results as needed (Needs citation)
def prokka(FASTAlist, FASTAfiles, PROKKAfiles, ORTHOfiles, CPUs):
    for file in FASTAlist:
        ProkkaPrefix = file.split(".")
        subprocess.Popen(["prokka", "--outdir", PROKKAfiles + file, "--compliant", "--prefix", ProkkaPrefix[0], "--force", "--genus", "Xanthomonus", "--species", "Oryzae", "--strain", ProkkaPrefix[0], "--cpus", CPUs, "--rfam", FASTAfiles + file], close_fds=True).communicate()[0]
        for item in os.listdir(PROKKAfiles + file + "/"):
            if item.endswith(".faa"):
                shutil.copyfile(PROKKAfiles + file + "/" + item, PROKKAfiles + "FAAs/" + item)
    os.system("grep '>' " + PROKKAfiles + "FAAs/*.faa | cut -d ' ' -f2- | sort | uniq -c > " + PROKKAfiles + "parsedProkka.txt")
    subprocess.Popen(["orthofinder", "-f", PROKKAfiles + "FAAs/", "-t", CPUs], close_fds=True).communicate()[0]
    for dir in os.listdir(PROKKAfiles + "FAAs/"):
        if "Results" in dir:
            for item in os.listdir(PROKKAfiles + "FAAs/" + dir):
                shutil.move(PROKKAfiles + "FAAs/" + dir + "/" + item, ORTHOfiles + item)


#Function for use of RVDminer, then passing of results to DisTAL
def RVDminer(RVDFiles, RVDsrc, RVDdst, DISdst):
    for file in RVDFiles:
        prefixRVD = file.split(".")
        subprocess.Popen(["perl", "RVDminer.pl", RVDsrc + file, prefixRVD[0]], close_fds=True).communicate()[0]
        for item in os.listdir(os.getcwd()):
            if prefixRVD[0] in item:
                shutil.move(os.getcwd() + "/" + item, RVDdst + item)
        os.system("cat " + RVDdst + "*TALS_aa.FASTA | sed 's/ Repeats//g' > " + DISdst + "rvdCombo.FASTA")
    shutil.move("./_Inline", RVDdst + "_Inline")
    comboFile = DISdst + "rvdCombo.FASTA"
    subprocess.Popen(["perl", "DisTAL_v1.3_Groups.pl", "-m", "T", comboFile, "disTALOut", "4.5"], close_fds=True).communicate()[0]
    shutil.move(os.getcwd() + "/" "Outputs", DISdst + "Outputs")


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
            else:
                print("Unable to find optimum value of K")


#Creates concatenated FAA file, adding file names to start of appropriate lines
def concatFaa(faaDir, RESULTSfiles):
    faaList = [file for file in os.listdir(faaDir) if file.lower().endswith(".faa")]
    faaList.sort()
    for file in faaList:
        prefix = file.split(".")
        with open(faaDir + file, 'r') as item:
            data = item.read()
        data = data.replace('>', '>' + prefix[0] + ' ')
        with open(faaDir + file, 'w') as item:
            item.write(data)
    with open(RESULTSfiles + "faaConcatenated.faa", 'wb') as outFile:
        for file in faaList:
            with open(faaDir + file, 'rb') as inFile:
                shutil.copyfileobj(inFile, outFile)
    return(RESULTSfiles + "faaConcatenated.faa")


#concatenates rvdNuc files to single CSV
def concatNuc(RVDfiles, RESULTSfiles):
    nucList = []
    for file in os.listdir(RVDfiles):
        if file.endswith(".csv"):
            nucList.append(RVDfiles + file)
    nucList.sort()
    combined = pd.concat([pd.read_csv(file, sep="\t") for file in nucList])
    combined.to_csv(RESULTSfiles + "rvdNucs.csv", index=False)
    return (RESULTSfiles + "rvdNucs.csv")


#Checks CSV files to be used for bound matrix for missing columns based on FASTAlist
def csvFix(Rfiles, FASTAlist):
    FASTAprefix = []
    for name in FASTAlist:
        temp = name.split(".")
        FASTAprefix.append(temp[0])
    if(os.path.isfile(Rfiles + "matTRF.csv")):
        matCSV = pd.read_csv(Rfiles + "matTRF.csv", index_col=0)
        matCols = list(matCSV.columns)
        matCompare = (list(set(FASTAprefix).difference(matCols)))
        for col in matCompare:
            print("Missing " + col + " in matTRF.csv. Adding, setting values to 0")
            matCSV[col] = int(0)
        matCSV.to_csv(Rfiles + "matTRF.csv")
    else:
        print("matTRF.csv was not generated, not applying csvFix")
    if(os.path.isfile(Rfiles + "OrthoGCMatrix.csv")):
        orthoCSV = pd.read_csv(Rfiles + "OrthoGCMatrix.csv", index_col=0)
        orthoCols = (list(orthoCSV.columns))
        orthoCompare = (list(set(FASTAprefix).difference(orthoCols)))
        for col in orthoCompare:
            print("Missing " + col + " in OrthoGCMatrix.csv. Adding, setting values to 0")
            orthoCSV[col] = int(0)
        orthoCSV.to_csv(Rfiles + "OrthoGCMatrix.csv")
    else:
        print("OrthoGCMatrix.csv was not generated, not applying csvFix")
    if(os.path.isfile(Rfiles + "distal_GroupMatrix.csv")):
        groupCSV = pd.read_csv(Rfiles + "distal_GroupMatrix.csv", index_col=0)
        groupCols = list(groupCSV.columns)
        groupCompare = (list(set(FASTAprefix).difference(groupCols)))
        for col in groupCompare:
            print("Missing " + col + " in distal_GroupMatrix.csv. Adding, setting values to 0")
            groupCSV[col] = int(0)
        groupCSV.to_csv(Rfiles + "distal_GroupMatrix.csv")
    else:
        print("distal_GroupMatrix.csv was not generated, not applying csvFix")


#Function for initializing scoary and second round of R modifications
def scoary(pipePath, providedCSV):
    shutil.copyfile(providedCSV, pipePath + "SCOARYfiles/providedCSV.csv")
    scorFile = (pipePath + "SCOARYfiles/providedCSV.csv")
    if scorFile.lower().endswith(".csv"):
        if(os.path.exists(pipePath + "Rfiles/boundMatrix.csv")):
            boundCSV = pd.read_csv(pipePath + "Rfiles/boundMatrix.csv", index_col=0)
            traitCSV = pd.read_csv(scorFile, index_col=0)
            boundCols = list(boundCSV.columns)
            traitRows = list(traitCSV.index)
            boundTraitCompare = list(set(boundCols) - set(traitRows))
            if 'V1' in boundTraitCompare:
                boundTraitCompare.remove('V1')
            if len(boundTraitCompare) is 0:
                return scorFile
            else:
                print("Columns of boundMatrix do not match rows of provided CSV, unable to run Scoary")
                return None
        else:
            print("Missing boundMatrix, unable to run Scoary")
            return None
    else:
        print("Provided csv does not end with csv extension, unable to proceed.")
        return None


#Function for routing needed data to R for parsing of various data sets
def scoaryParse(SCOARYfiles, Rfiles, DISTALfiles, TRFfiles, ORTHOfiles, RESULTSfiles, Logging):
    for file in os.listdir(SCOARYfiles):
        if file.endswith("results.csv"):
            scoaryCSV = (SCOARYfiles + file)
    try:
        scoaryCSV
    except:
        print("scoary CSV doesn't exist, unable to continue process")
        sys.exit()
    else:
        repeatsCSV = Rfiles + "RepeatNames.csv"
        distalGroups = DISTALfiles + "Outputs/disTALOut.TALgroups.csv"
        parsedTRF = TRFfiles + "trfParsed.txt"
        orthogroupsTXT = ORTHOfiles + "Orthogroups.txt"
        comboFile = DISTALfiles + "rvdCombo.FASTA"
        faaFile = (RESULTSfiles + "faaConcatenated.faa")
        rvdNucs = (RESULTSfiles + "rvdNucs.csv")
        subprocess.Popen(["Rscript", "addScripts/BactRThree.R", scoaryCSV, repeatsCSV, distalGroups, parsedTRF, orthogroupsTXT, comboFile, RESULTSfiles, faaFile, rvdNucs, Logging], close_fds=True).communicate()[0]


#Parsing of BayesGenerator.csv, then shuttling of data to Bayes
def bayesPool(pipePath):
    bayesGen = (pipePath + "BAYESfiles/BayesGenerator.csv")
    if(os.path.exists(bayesGen)):
        if(os.path.exists(pipePath + "BAYESfiles/Ktree.nexus")):
            bayesData = pd.read_csv(bayesGen, index_col=0)
            with open(pipePath + "Results/bayesResults.txt", 'w') as bayesOut:
                bayesOut.writelines("Gene\tDependent\tIndependent\t2*(dep-ind)\n")
            bayesList = (bayesData.columns)
            pool = mp.Pool(processes=10,)
            pool.starmap(bayesCall, zip(bayesList, repeat(pipePath)))
            pool.close()
            pool.join()
        else:
            print("Missing Ktree.nexus file, unable to process Bayes.")
    else:
        print("Missing BayesGenerator.csv, unable to process Bayes.")

#Helper function for pool process of BAYESfiles
def bayesCall(iterableColumns, path):
    nexus = (path + "BAYESfiles/Ktree.nexus")
    bayesData = pd.read_csv(path + "BAYESfiles/BayesGenerator.csv", index_col=0)
    traitData = pd.read_csv(path + "SCOARYfiles/providedCSV.csv", index_col=0)
    df = bayesData[[iterableColumns]].join(traitData)
    dfName = (path + "BAYESfiles/" + iterableColumns + "_mat.txt")
    df.to_csv(dfName, header=False, sep= "\t")

    #Subprocess call for dependent Bayes
    depProc = subprocess.Popen(["BayesTraitsV3", nexus, dfName], stdin=subprocess.PIPE, text=True)
    depCom=("3\n2\nPriorAll exp 10\nStones 100 500\nCor 1\nRun\n")
    depProc.stdin.write(depCom)
    depProc.stdin.flush()
    depProc.communicate()[0]

    if(os.path.isfile(dfName + ".Stones.txt")):
        os.rename(dfName + ".Stones.txt", path + "BAYESfiles/" + iterableColumns + ".dependant")
        with open(path + "BAYESfiles/" + iterableColumns + ".dependant", 'r') as depFile:
            content = depFile.readlines()
            for line in content:
                if "Log marginal likelihood" in line:
                    temp = line.split()
                    tempDep = temp[3]
    else:
        tempDep = 'NA'


    #Subprocess call for independent Bayes
    indProc = subprocess.Popen(["BayesTraitsV3", nexus, dfName], stdin=subprocess.PIPE, text=True)
    indCom=("2\n2\nPriorAll exp 10\nStones 100 500\nCor 1\nRun\n")
    indProc.stdin.write(indCom)
    indProc.stdin.flush()
    indProc.communicate()[0]

    if(os.path.isfile(dfName + ".Stones.txt")):
        os.rename(dfName+ ".Stones.txt", path + "BAYESfiles/" + iterableColumns + ".independant")
        with open(path + "BAYESfiles/" + iterableColumns + ".independant", 'r') as indepFile:
            content = indepFile.readlines()
            for line in content:
                if "Log marginal likelihood" in line:
                    temp = line.split()
                    tempInd = temp[3]
    else:
        tempInd = "NA"

    #Check variables for dep and ind
    try:
        tempDep
    except NameError:
        print("Missing value for " + iterableColumns + " in dependant file.")
        dependantVal = 'NA'
    else:
        dependantVal = tempDep
    try:
        tempInd
    except NameError:
        print("Missing value for " + iterableColumns + " in independant file.")
        independantVal = 'NA'
    else:
        independantVal = tempInd

    #Get combo value, 2*(dependant - independant)
    if(dependantVal != "NA" and independantVal != 'NA'):
        finalVal = (2*(float(dependantVal) - float(independantVal)))
    else:
        finalVal = 'NA'

    #Write variables to results file
    with open(path + "Results/bayesResults.txt", 'a') as bayesOut:
        bayesOut.writelines(iterableColumns + "\t" + dependantVal + "\t" + independantVal + "\t" + str(finalVal) + "\n")
    for item in os.listdir(path + "BAYESfiles/"):
        if iterableColumns in item:
            os.remove(path + "BAYESfiles/" + item)
