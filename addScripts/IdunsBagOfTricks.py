import subprocess
import sys
import shutil
import os
import re

sys.dont_write_bytecode = True



#Run FASTA files through prokka, parse FAA'srun results through Ortho, move results as needed (Needs citation)
def prokka(files, src, dst, ortho, cpus):
    for file in files:
        prefix = file.split(".")
        with open(dst + 'prokkaErr.txt', 'w') as prokkaErr:
            subprocess.Popen(["prokka", "--outdir", dst + file, "--compliant", "--prefix", prefix[0], "--force", "--genus", "Xanthomonus", "--species", "Oryzae", "--strain", prefix[0], "--cpus", "4", "--rfam", src + file], stderr=prokkaErr, close_fds=True).communicate()[0]
        for item in os.listdir(dst + file + "/"):
            if item.endswith(".faa"):
                shutil.copyfile(dst + file + "/" + item, dst + "FAAs/" + item)
    os.system("grep '>' " + dst + "FAAs/*.faa | cut -d ' ' -f2- | sort | uniq -c > " + dst + "parsedProkka.txt")
    with open(ortho + "orthoErr.txt", 'w') as orthoErr:
        subprocess.Popen(["orthofinder", "-f", dst + "FAAs/", "-t", cpus], close_fds=True, stderr=orthoErr).communicate()[0]
    for dir in os.listdir(dst + "FAAs/"):
        if "Results" in dir:
            for item in os.listdir(dst + "FAAs/" + dir):
                shutil.move(dst + "FAAs/" + dir + "/" + item, ortho + item)



#Function for use of RVDminer, then passing of results to DisTAL
def RVDminer(files, src, RVDdst, DISdst):
    for file in files:
        prefix = file.split(".")
        with open(RVDdst + "rvdErr.txt", 'w') as rvdErr:
            subprocess.Popen(["perl", "RVDminer.pl", src + file, prefix[0]], close_fds=True, stderr=rvdErr).communicate()[0]
        for item in os.listdir(os.getcwd()):
            if prefix[0] in item:
                shutil.move(os.getcwd() + "/" + item, RVDdst + item)
        os.system("cat " + RVDdst + "*TALS_aa.FASTA | sed 's/ Repeats//g' > " + DISdst + "rvdCombo.FASTA")
    shutil.move("./_Inline", RVDdst + "_Inline")
    comboFile = DISdst + "rvdCombo.FASTA"
    with open(DISdst + 'disErr', 'w') as disErr:
        subprocess.Popen(["perl", "DisTAL_v1.3_Groups.pl", "-m", "T", comboFile, "disTALOut", "4.5"], close_fds=True, stderr=disErr).communicate()[0]
    shutil.move(os.getcwd() + "/" "Outputs", DISdst + "Outputs")



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


#Parsing of TRF files to single text file
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



#Function for use of KSNP3, appropriate parsing, and movement of files
def ksnpCall(faPath, ksnpPath, faList, CPUs):
    ksnpGenomes = ksnpPath + "ksnpGenomes.txt"
    with open(ksnpGenomes, "w") as outFile:
        for file in faList:
            prefix = file.split(".")
            outFile.write(os.path.abspath(faPath + file) + "\t" + prefix[0] + "\n")
    subprocess.Popen(["MakeFasta", ksnpGenomes, ksnpPath + "ForKchoser"], close_fds=True).communicate()[0]
    subprocess.Popen(["Kchooser", ksnpPath + "ForKchoser"], close_fds=True).communicate()[0]
    shutil.move("Kchooser.report", ksnpPath + "Kchooser.report")
    with open(ksnpPath + "Kchooser.report", "r") as kOut:
        content = kOut.readlines()
        for line in content:
            if re.match(r'The optimum value of K is (.*)\.', line):
                tempInt = re.findall(r'The optimum value of K is (.*)\.', line)
                ksnpInt = tempInt[0]
    subprocess.Popen(["kSNP3", "-in", ksnpGenomes, "-k", str(ksnpInt), "-outdir", ksnpPath + "kSNP3_results", "-ML", "-CPU", str(CPUs)], close_fds=True).communicate()[0]
    shutil.move("fasta_list", ksnpPath + "fasta_list")
