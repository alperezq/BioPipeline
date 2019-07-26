#!/usr/bin/python
#Bioinformatics Pipeline
from addScripts import IdunsProkkaOrtho as IPO
from addScripts import IdunsRVDtoDIS as IRD
from addScripts import IdunsKSNP3 as IK3
import os, shutil #Necessary for directory checks and file movements
import time #To breakup processes and make output readable
import sys #Allows exiting of code in case of irreconcilable error
import argparse #Allows use of command line style argument call
import subprocess #For execution of outside scripts
import multiprocessing as mp #Run processes & run multiple processes at once
import pandas as pd
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

#Main section / execution of code

#Pipe variables
projName = args.name
processors = str(args.processors)

#Processor variable, halves if greater than 10 per project guidelines by Alvaro
if processors <= 10:
    CPUs = processors
else:
    CPUs = (int(processors)/2)
CPUs = str(CPUs)

#Initiate path variables and file location variables (created after start of main section)
pipePath, fastaPath = pipeStart(args.name, args.fasta)
FASTAfiles = pipePath + "FASTAfiles/"
TRFfiles = pipePath + "TRFfiles/"
PROKKAfiles = pipePath + "PROKKAfiles/"
ORTHOfiles = pipePath + "ORTHOfiles/"
RVDfiles = pipePath + "RVDfiles/"
DISTALfiles = pipePath + "DISTALfiles/"
Rfiles = pipePath + "Rfiles/"
KSNP3files = pipePath + "KSNP3files/"
SCOARYfiles = pipePath + "SCOARYfiles/"
RESULTSfiles = pipePath + "Results/"
repeatCSV = Rfiles + "RepeatNames.csv"
boundCSV = Rfiles + "boundMatrix.csv"
providedCSV = args.scoary
talGroupsCSV = DISTALfiles + "Outputs/disTALout.TALgroups.csv"
trfTXT = TRFfiles + "trfParsed.txt"
orthogroupsTXT = ORTHOfiles + "Orthogroups.txt"
kTree = KSNP3files + "kSNP3_results/tree.ML.tre"
rvdFASTA = DISTALfiles + "rvdCombo.FASTA"

#Gather FASTA files, copy to project folder for further use
genomeNumber, FASTAlist = collectFasta(fastaPath, FASTAfiles)

#Call to TandemRepeatsFinder, done individually to allow processing of large batches of Files
if __name__== '__main__':
    workerPool = mp.Pool(processes=int(CPUs))
    workerPool.map(tandemRepeatFinder, FASTAlist)
    workerPool.close()
IPO.trfParse(TRFfiles, FASTAlist)

#Establish first set of processes for the pipeline and pass their relevant parameters
prokkaProcess = mp.Process(target = IPO.prokka, args =(FASTAlist, FASTAfiles, PROKKAfiles, ORTHOfiles, CPUs,))
RVDProcess = mp.Process(target = IRD.RVDminer, args = (FASTAlist, FASTAfiles, RVDfiles, DISTALfiles,))
KSNP3Process = mp.Process(target = IK3.ksnpCall, args = (FASTAfiles, KSNP3files, FASTAlist, CPUs,))

#Start processes
prokkaProcess.start()
RVDProcess.start()
KSNP3Process.start()

#Rejoin processes with main thread, won't continue till each finishes
prokkaProcess.join()
RVDProcess.join()
KSNP3Process.join()

#Call R script for further parsing of data
subprocess.Popen(["Rscript", "addScripts/IdunsRScript.R", pipePath], close_fds=True).communicate()[0]

#Call Scoary if it is supplied the necessary CSV, compares rows of CSV with colums of boundMatrix.csv first
if providedCSV is not None:
    if args.scoary.lower().endswith(".csv"):
        boundCSV = pd.read_csv(Rfiles + "boundMatrix.csv", index_col=0)
        traitCSV = pd.read_csv(args.scoary, index_col=0)
        boundCols = list(boundCSV.columns)
        traitRows = list(traitCSV.index)
        boundTraitCompare = len(list(set(boundCols) - set(traitRows)))
        if boundTraitCompare is 0:
            subprocess.Popen(["scoary", "-t", providedCSV, "-g", Rfiles + "boundMatrix.csv", "-s", "2", "-o", SCOARYfiles], close_fds=True).communicate()[0]
            IK3.ksnpParse(SCOARYfiles, RFiles, providedCSV, DISTALfiles, TRFfiles, ORTHOfiles, KSNP3files, RESULTSfiles)
        else:
            print("Columns of boundMatrix do not match rows of provided CSV, unable to run Scoary")
