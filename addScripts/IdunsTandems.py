import subprocess
import sys
import shutil
import os


sys.dont_write_bytecode = True


#Run FASTA files through Tandem Repeat Finder
#Citation:          G. Benson,
#                   "Tandem repeats finder: a program to analyze DNA sequences"
#                   Nucleic Acids Research (1999)
#                   Vol. 27, No. 2, pp. 573-580.
def tandemRepeatFinder(TandemSrc, TandemDst, TandemFiles):
    for initFile in TandemFiles:
        subprocess.Popen(["TandemRepeatsFinder", TandemSrc + initFile, "2", "7", "7", "80", "10", "50", "500", "-f", "-h"], close_fds=True).communicate()[0]
    for datFile in os.listdir(os.getcwd()):
        print datFile
        for endFile in TandemFiles:
            if endFile in datFile:
                shutil.move(os.getcwd() + "/" + datFile, TandemDst + datFile)
    trfParse(TandemDst)


#Parsing of TRF files to single text file
def trfParse(TanParseDir):
    with open(TanParseDir + "trfParsed.txt", "w") as outFile:
        for file in os.listdir(TanParseDir):
            if file.endswith(".dat"):
                prefixTRF = file.split(".")[0]
                with open(TanParseDir + file) as inFile:
                    for line in inFile.readlines():
                        if (len(line) >= 45) and ("Sequence" not in line):
                            entries = line.split()
                            entry = (prefixTRF + " " + entries[0] + " " + entries[1] + " " + entries[2] + " " + entries[3] + " " + entries[4] + " " + entries[13] + " " + entries[14] + "\n")
                            outFile.write(entry)
