import subprocess

#Run FASTA files through prokka, set to run with 4 CPUS, may add ability to modify at a later date (Needs citation)
def prokka(files, src, dst):
    for file in files:
        prefix = file.split(".")
        subprocess.Popen(["prokka", "--outdir", dst + file, "--compliant", "--prefix", prefix[0], "--force", "--genus", "Xanthomonus", "--species", "Oryzae", "--strain", prefix, "--cpus", "4", "--rfam", src + file], close_fds=True).communicate()[0]
    subprocess.Popen(["orthofinder", "-f", dst + "*/*.faa"], close_fds=True).communicate()[0]
