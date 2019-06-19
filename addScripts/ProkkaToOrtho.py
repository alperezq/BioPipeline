import subprocess
import sys
import shutil
import os

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
