import subprocess
import sys
import shutil
import os


sys.dont_write_bytecode = True


#Run FASTA files through prokka, parse FAA'srun results through Ortho, move results as needed (Needs citation)
def prokka(ProkkaFiles, ProkkaSrc, ProkkaDst, ortho, ProkkaOrthoCpus):
    for file in ProkkaFiles:
        ProkkaPrefix = file.split(".")
        subprocess.Popen(["prokka", "--outdir", ProkkaDst + file, "--compliant", "--prefix", ProkkaPrefix[0], "--force", "--genus", "Xanthomonus", "--species", "Oryzae", "--strain", prefix[0], "--cpus", ProkkaCpus, "--rfam", ProkkaSrc + file], close_fds=True).communicate()[0]
        for item in os.listdir(ProkkaDst + file + "/"):
            if item.endswith(".faa"):
                shutil.copyfile(ProkkaDst + file + "/" + item, ProkkaDst + "FAAs/" + item)
    os.system("grep '>' " + ProkkaDst + "FAAs/*.faa | cut -d ' ' -f2- | sort | uniq -c > " + ProkkaDst + "parsedProkka.txt")
    subprocess.Popen(["orthofinder", "-f", ProkkaDst + "FAAs/", "-t", ProkkaOrthoCpus], close_fds=True).communicate()[0]
    for dir in os.listdir(ProkkaDst + "FAAs/"):
        if "Results" in dir:
            for item in os.listdir(ProkkaDst + "FAAs/" + dir):
                shutil.move(ProkkaDst + "FAAs/" + dir + "/" + item, ortho + item)
