#!/usr/bin/python3.7
#Bioinformatics Pipeline
from addScripts import IdunsDawn as IDawn
from addScripts import IdunsProkkaOrtho as IPO
from addScripts import IdunsRVDtoDIS as IRD
from addScripts import IdunsKSNP3 as IK3
from addScripts import IdunsBridge as IBridge
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
parserPS.add_argument("name", help="Name of project. Must be non-existant directory if running first half of pipeline. If just running the second half, must use project directory created by / in same format as Idun's Pipeline")
parserPS.add_argument("processors", help="Number of processors to utilize for this job", type=int)
parserPS.add_argument("-P","--pipeStart", nargs="?", help="Optional argument to run first half of pipeline. Must provide directory with FASTA files for processing")
parserPS.add_argument("-S","--scoary", nargs='?', help="Optional argument to add a CSV for Scoary processing. Scoary, and second portion of pipeline, will not run without this, and columns must match boundmatrix.csv that is generated")
args = parserPS.parse_args()

def tandemRepeatFinder(tanFile):
    subprocess.Popen(["TandemRepeatsFinder", FASTAfiles + tanFile, "2", "7", "7", "80", "10", "50", "500", "-f", "-h"], close_fds=True).communicate()[0]

#Main section / execution of code
if __name__== '__main__':

    #Assignment of pipe scope variables
    pipeStart = False
    pipeScoary = False
    if args.scoary is not None:
        pipeScoary = True
    if args.pipeStart is not None:
        pipeStart = True

    #Processor variable, halves if greater than 10 per project guidelines by Alvaro
    if int(args.processors) <= 10:
        CPUs = args.processors
    else:
        CPUs = (int(args.processors)//2)
    CPUs = str(CPUs)

    #If neither pipeStart or pipeScoary is initiated with True, pipeline won't run
    if pipeStart is False and pipeScoary is False:
        print("Neither start or end of pipe initialized. Exiting...")
        sys.exit()

    #Assignment of pipePath and FASTAlist, dependent on mode being run
    if pipeStart is True:
        if args.pipeStart is None:
            print("Must provide a directory with FASTA files when running initial section of pipeline.\nExiting...")
            sys.exit()
        pipePath, fastaPath = IDawn.pipeStart(args.name, args.pipeStart)
    else:
        pipePath, FASTAlist = IDawn.pipeDetour(args.name)

    #Initiate path variables and file location variables (created after start of main section)
    FASTAfiles = pipePath + "FASTAfiles/"
    TRFfiles = pipePath + "TRFfiles/"
    PROKKAfiles = pipePath + "PROKKAfiles/"
    ORTHOfiles = pipePath + "ORTHOfiles/"
    RVDfiles = pipePath + "RVDfiles/"
    DISTALfiles = pipePath + "DISTALfiles/"
    Rfiles = pipePath + "Rfiles/"
    KSNP3files = pipePath + "KSNP3files/"
    SCOARYfiles = pipePath + "SCOARYfiles/"
    BAYESfiles = pipePath + "BAYESfiles/"
    RESULTSfiles = pipePath + "Results/"
    LOGfiles = pipePath + "Logging/"
    providedCSV = args.scoary

    #Gather FASTA files, copy to project folder for further use
    if pipeStart is True:
        FASTAlist = IDawn.collectFasta(fastaPath, FASTAfiles)

        #Call to TandemRepeatsFinder, done individually to allow processing of large batches of Files
        workerPool = mp.Pool(processes=int(CPUs),)
        workerPool.map(tandemRepeatFinder, FASTAlist)
        workerPool.close()
        workerPool.join()

        #Parsing of TRF files
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

        #Creation of faaConcatenated and rvdNucs files for end Results
        IPO.concatFaa(PROKKAfiles + "FAAs/", RESULTSfiles)
        IRD.concatNuc(RVDfiles, RESULTSfiles)

    #Second section of pipeline, requires secondHalf to be True
    if pipeScoary is True:

        #Call R script for further parsing of data
        subprocess.Popen(["Rscript", "addScripts/IdunsRScript.R", pipePath], close_fds=True).communicate()[0]
        IBridge.csvFix(Rfiles, FASTAlist)
        subprocess.Popen(["Rscript", "addScripts/IdunsRBridge.R", pipePath], close_fds=True).communicate()[0]

        #Call Scoary if it is supplied the necessary CSV, compares rows of CSV with colums of boundMatrix.csv first
        if providedCSV is not None:
            scorFile = IBridge.scoary(pipePath, providedCSV)
            subprocess.Popen(["Rscript", "addScripts/Iduns3rdR.R", pipePath], close_fds=True).communicate()[0]
            subprocess.Popen(["scoary", "-t", scorFile, "-g", pipePath + "Rfiles/boundMatrix.csv", "-s", "2", "-o", pipePath + "SCOARYfiles/"], close_fds=True).communicate()[0]
            IK3.ksnpParse(SCOARYfiles, Rfiles, scorFile, DISTALfiles, TRFfiles, ORTHOfiles, KSNP3files, RESULTSfiles, RVDfiles, PROKKAfiles + "FAAs/")

            #Call BayesTraitsV3 on prior results of pipeline
            #IBridge.bayesPool(pipePath)
