import subprocess

def TEST():
    subprocess.Popen(["perl", "RVDminer/RVDminer.pl"], close_fds=True).communicate()[0]

print("Attempting Test")
TEST()
print("Attempted Test")
