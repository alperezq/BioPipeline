import subprocess
import sys
import shutil
import os


sys.dont_write_bytecode = True


#Run FASTA files through prokka, parse FAA'srun results through Ortho, move results as needed (Needs citation)
def prokka(ProkkaFiles, ProkkaSrc, ProkkaDst, ortho, ProkkaOrthoCpus):
    for file in ProkkaFiles:
        ProkkaPrefix = file.split(".")
        subprocess.Popen(["prokka", "--outdir", ProkkaDst + file, "--compliant", "--prefix", ProkkaPrefix[0], "--force", "--genus", "Xanthomonus", "--species", "Oryzae", "--strain", ProkkaPrefix[0], "--cpus", ProkkaOrthoCpus, "--rfam", ProkkaSrc + file], close_fds=True).communicate()[0]
        for item in os.listdir(ProkkaDst + file + "/"):
            if item.endswith(".faa"):
                shutil.copyfile(ProkkaDst + file + "/" + item, ProkkaDst + "FAAs/" + item)
    os.system("grep '>' " + ProkkaDst + "FAAs/*.faa | cut -d ' ' -f2- | sort | uniq -c > " + ProkkaDst + "parsedProkka.txt")
    subprocess.Popen(["orthofinder", "-f", ProkkaDst + "FAAs/", "-t", ProkkaOrthoCpus], close_fds=True).communicate()[0]
    for dir in os.listdir(ProkkaDst + "FAAs/"):
        if "Results" in dir:
            for item in os.listdir(ProkkaDst + "FAAs/" + dir):
                shutil.move(ProkkaDst + "FAAs/" + dir + "/" + item, ortho + item)


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
