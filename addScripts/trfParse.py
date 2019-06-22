import os
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True

def trfParse(dir):
    with open(dir + "trfParsed.txt", "w") as outFile:
        for file in os.listdir(dir):
            if file.endswith(".dat"):
                prefix = file.split(".")[0]
                with open(dir + file) as inFile:
                    for line in inFile.readlines():
                        if (len(line) >= 45) and ("Sequence" not in line):
                            entries = line.split()
                            entry = (prefix + " " + entries[0] + " " + entries[1] + " " + entries[2] + " " + entries[3] + " " + entries[4] + " " + entries[13] + " " + entries[14] + "\n")
                            outFile.write(entry)

#Run FASTA files through Tandem Repeat Finder
#Citation:          G. Benson,
#                   "Tandem repeats finder: a program to analyze DNA sequences"
#                   Nucleic Acids Research (1999)
#                   Vol. 27, No. 2, pp. 573-580.
def tandemRepeatFinder(src, dst, files):
    for initFile in files:
        subprocess.Popen(["TandemRepeatsFinder", src + initFile, "2", "7", "7", "80", "10", "50", "500", "-f", "-h"], close_fds=True).communicate()[0]
    for datFile in os.listdir(os.getcwd()):
        print datFile
        for endFile in files:
            print endFile
            if endFile in datFile:
                print(datFile)
                shutil.move(os.getcwd() + "/" + datFile, dst + datFile)
    trfParse(dst)
