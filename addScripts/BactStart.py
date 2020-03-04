#!/usr/bin/env python3
import os
import sys
import subprocess
import shutil

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
            os.mkdir(fullPath + "BAYESfiles")
            os.mkdir(fullPath + "Results")
            os.mkdir(fullPath + "Logging")
        except OSError:
            print("Creation of the directory failed. Exiting...")
            sys.exit()
        else:
            print("Succesfully created the directory")
    return fullPath, filePath;

def pipeDetour(name):
    if not name.endswith("/"):
        name = name + "/"
    fileExt = [".fa", ".fasta", ".fas", "fna"]
    fileList = []
    for file in os.listdir(name + "FASTAfiles/"):
        if file.lower().endswith(tuple(fileExt)):
            fileList.append(file)
    return name, fileList;


#Checks FASTA path for .FASTA txt files, creates new directory in project directory with copies of files
def collectFasta(src, dst):
    numberOfFiles = 0
    fileExt = [".fa", ".fasta", ".fas", "fna"]
    fileList = []
    for file in os.listdir(src):
        if file.lower().endswith(tuple(fileExt)):
            fileList.append(file)
            shutil.copyfile(src + file, dst + file)
    return fileList;
