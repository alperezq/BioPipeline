#!/usr/bin/env python3
import subprocess
import sys
import shutil
import os


sys.dont_write_bytecode = True


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
