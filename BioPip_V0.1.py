#!/usr/bin/python
#Bioinformatics Pipeline
import ProkkaToOrtho as PTO #Import Py script for Prokka and Orthofinder call
import RVDtoDIS as RtoD #Import Py script for RVDminer and DisTAL call
import os, shutil #Necessary for directory checks and file movements
import time #To breakup processes and make output readable
import sys #Allows exiting of code in case of irreconcilable error
import glob #Manipulation of files in directories
import argparse #Allows use of command line style argument call
import subprocess #For execution of outside scripts
from multiprocessing import Pool #For concurrent processes

#Parser from argparse for command line refinement
parserPS = argparse.ArgumentParser()
parserPS.add_argument("name", help="Name of project, must be non-existant directory")
parserPS.add_argument("fasta", help="Full valid directory path where FASTA files are stored. Files must be text files with .FASTA extension")
parserPS.add_argument("processors", help="Number of processors to utilize for this job", type=int)
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
            os.mkdir(fullPath + "ORTHOfiles")
            os.mkdir(fullPath + "RVDfiles")
            os.mkdir(fullPath + "DISTALfiles")
        except OSError:
            print("Creation of the directory failed. Exiting...")
            sys.exit()
        else:
            print("Succesfully created the directory")
    return fullPath, filePath;

#Checks FASTA path for .FASTA txt files, creates new directory in project directory with copies of files
def collectFasta(src, dst):
    numberOfFiles = 0
    fileList = []
    for file in os.listdir(src):
        if file.endswith(".fa"):
            numberOfFiles+=1
            fileList.append(file)
            shutil.copyfile(src + "/" + file, dst + "/" + file)
    return numberOfFiles, fileList;

#Run FASTA files through Tandem Repeat Finder
#Citation:          G. Benson,
#                   "Tandem repeats finder: a program to analyze DNA sequences"
#                   Nucleic Acids Research (1999)
#                   Vol. 27, No. 2, pp. 573-580.
def tandemRepeatFinder(files):
    for initFile in files:
        subprocess.Popen(["trf", fastaPath + initFile, "2", "7", "7", "80", "10", "50", "500", "-f", "-h"], close_fds=True).communicate()[0]
    for datFile in os.listdir(os.getcwd()):
        for endFile in files:
            if endFile in datFile:
                shutil.move(os.getcwd() + "/" + datFile, TRFfiles + datFile)

#Pipe variables
projName = args.name
processors = args.processors

#Main section / execution of code
#Initiate path variables
pipePath, fastaPath = pipeStart(args.name, args.fasta)
FASTAfiles = pipePath + "FASTAfiles/"
TRFfiles = pipePath + "TRFfiles/"
PROKKAfiles = pipePath + "PROKKAfiles/"
ORTHOfiles = pipePath + "ORTHOfiles/"
RVDfiles = pipePath + "RVDfiles/"
DISTALfiles = pipePath + "DISTALfiles/"

#Gather FASTA files, copy to project folder for further use
genomeNumber, FASTAlist = collectFasta(fastaPath, FASTAfiles)

#Run Tandem Repeat Finder on FASTA files
tandemRepeatFinder(FASTAlist)

#Call to ProkkaToOrtho script
PTO.prokka(FASTAlist, FASTAfiles, PROKKAfiles)

#Call to RVDtoDIS script
RtoD.RVDminer(FASTAlist, FASTAfiles, RVDfiles, DISTALfiles)
