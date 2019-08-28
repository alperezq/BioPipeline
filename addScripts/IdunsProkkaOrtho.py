#!/usr/bin/env python3
import subprocess
import sys
import shutil
import os


sys.dont_write_bytecode = True


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
