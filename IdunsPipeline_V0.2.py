#!/usr/bin/python
#Bioinformatics Pipeline
from addScripts import IdunsBagOfTricks as IDBT
import os, shutil #Necessary for directory checks and file movements
import time #To breakup processes and make output readable
import sys #Allows exiting of code in case of irreconcilable error
import argparse #Allows use of command line style argument call
import subprocess #For execution of outside scripts
import multiprocessing as mp #Run proceses & run multiple processes at once
import time

#Prevents creation of .pyc files when running
sys.dont_write_bytecode = True

#Parser from argparse for command line refinement
parserPS = argparse.ArgumentParser()
parserPS.add_argument("name", help="Name of project, must be non-existant directory")
parserPS.add_argument("fasta", help="Full valid directory path where FASTA files are stored. Files must be text files with .FASTA extension")
parserPS.add_argument("processors", help="Number of processors to utilize for this job", type=int)
parserPS.add_argument("-s","--scoary", nargs='?', help="Optional argument to add a CSV for Scoary processing. Scoary will not run without this, and columns must match boundmatrix.csv that is generated")
args = parserPS.parse_args()

#Verify + Create directories
def pipeStart(name, filePath):
    if not name.endswith("/"):
        name = name + "/"
    fullPath = "./" + name
    if not filePath.endswith("/"):
        filePath = filePath + "/"
    if not os.path.isdir(filePath):
        print ("%s is not a valid directory" % filePath)
        sys.exit()
    if os.path.isdir(fullPath):
        print ("The directory %s already exists." % fullPath)
        sys.exit()
    else:
        print ("Creating directory %s" % fullPath)
        try:
            os.mkdir(fullPath)
            os.mkdir(fullPath + "FASTAfiles")
            os.mkdir(fullPath + "TRFfiles")
            os.mkdir(fullPath + "PROKKAfiles")
            os.mkdir(fullPath + "PROKKAfiles/FAAs")
            os.mkdir(fullPath + "ORTHOfiles")
            os.mkdir(fullPath + "RVDfiles")
            os.mkdir(fullPath + "DISTALfiles")
            os.mkdir(fullPath + "Rfiles")
            os.mkdir(fullPath + "KSNP3files")
            os.mkdir(fullPath + "SCOARYfiles")
            os.mkdir(fullPath + "Results")
        except OSError:
            print("Creation of the directory failed. Exiting...")
            sys.exit()
        else:
            print("Succesfully created the directory")
    return fullPath, filePath;

#Checks FASTA path for .FASTA txt files, creates new directory in project directory with copies of files
def collectFasta(src, dst):
    numberOfFiles = 0
    fileExt = [".fa", ".fasta", ".fas", "fna"]
    fileList = []
    for file in os.listdir(src):
        if file.lower().endswith(tuple(fileExt)):
            numberOfFiles+=1
            fileList.append(file)
            shutil.copyfile(src + file, dst + file)
    return numberOfFiles, fileList;

#Pipe variables
projName = args.name
processors = str(args.processors)

if processors <= 10:
    CPUs = 10
else:
    CPUs = (processors/2)

#Main section / execution of code
#Initiate path variables
pipePath, fastaPath = pipeStart(args.name, args.fasta)
FASTAfiles = pipePath + "FASTAfiles/"
TRFfiles = pipePath + "TRFfiles/"
PROKKAfiles = pipePath + "PROKKAfiles/"
ORTHOfiles = pipePath + "ORTHOfiles/"
RVDfiles = pipePath + "RVDfiles/"
DISTALfiles = pipePath + "DISTALfiles/"
KSNP3files = pipePath + "KSNP3files/"
SCOARYfiles = pipePath + "SCOARYfiles/"
RESULTSfiles = pipepath + "Results/"

#Gather FASTA files, copy to project folder for further use
genomeNumber, FASTAlist = collectFasta(fastaPath, FASTAfiles)

#Establish first set of processes for the pipeline and pass their relevant parameters
tandemProcess = mp.Process(target = IDBT.tandemRepeatFinder, args = (fastaPath, TRFfiles, FASTAlist,))
prokkaProcess = mp.Process(target = IDBT.prokka, args =(FASTAlist, FASTAfiles, PROKKAfiles, ORTHOfiles, CPUs))
RVDProcess = mp.Process(target = IDBT.RVDminer, args = (FASTAlist, FASTAfiles, RVDfiles, DISTALfiles,))
KSPN3Process = mp.process(target = IDBT.ksnpCall, args = (FASTAfiles, KSNP3files, FASTAlist, CPUs))

#Start processes
tandemProcess.start()
prokkaProcess.start()
RVDProcess.start()
KSPN3Process.start()

#Rejoin processes with main thread, won't continue till each finishes
tandemProcess.join()
prokkaProcess.join()
RVDProcess.join()
KSNP3Process.join()

#Call R script for further parsing of data
subprocess.Popen(["Rscript", "addScripts/IdunsRScript.R", pipePath], close_fds=True).communicate()[0]

#Call Scoary if it is supplied the necessary CSV
#Add check if the file is a CSV or not, then copy CSV to SCOARYfiles, add tree file
if args.scoary is not None:
    if arg.scoary.endswith(.csv):
        subprocess.Popen(["Scoary", "-t", args.scoary, "-g", Rfiles + "boundmatrix.csv", "-s", "-2", "-o", SCOARYfiles], close_fds=True).communicate()[0]
