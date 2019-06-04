import subprocess
import os
import shutil


def RVDminer(files, src, dst, disOUT):
    for file in files:
        prefix = file.split(".")
        subprocess.Popen(["perl", "RVDminer/RVDminer.pl", src + file, prefix], close_fds=True).communicate()[0]
    for item in os.listdir(os.getcwd):
        for RVOUTfile in files:
            if RVOUTfile in item:
                shutil.move(os.getcwd() + "/" + item, dst + item)
    subprocess.Popen(["cat " + dst + "*TALS_aa.FASTA > " + dst + "rvdOut.txt"], close_fds=True).communicate()[0]
    comboFile = dst + "rvdOut.txt"
    subprocess.Popen(["perl", "DisTAL/DisTAL_v1.3_Groups.pl", "-m T", comboFile, disOUT, "4.5"], close_fds=True).communicate()[0]
    
