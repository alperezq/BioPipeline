import subprocess
import os
import shutil
import sys

sys.dont_write_bytecode = True

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
    shutil.move(os.getcwd() + "/Outputs", DISdst + "Outputs")
