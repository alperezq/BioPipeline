#!/usr/bin/env Rscript
#Declaration of libraries
library(ggplot2)
library(reshape2)
library(ape)
library(stringr)
library(plyr)
library(dplyr)
library(argparse)

#Take in args from IdunsBridge
args <- commandArgs(trailingOnly=TRUE)
projectPath = args[1]

#Assign paths to variables
orthoGCMatrix <- read.csv(paste(projectPath, "Rfiles/OrthoGCMatrix.csv", sep=""))
matTRF <- read.csv(paste(projectPath, "Rfiles/matTRF.csv", sep=""))
groupMat <- read.csv(paste(projectPath, "Rfiles/distal_GroupMatrix.csv", sep=""))

#Creation of bound matrix after manipulations in python3
boundMatrix <- bind_rows(orthoGCMatrix, matTRF, groupMat)
boundMatrixFix <- slice(boundMatrix, 1:(n()-2))
colnames(boundMatrixFix) <- gsub(".","-", colnames(boundMatrixFix), fixed = TRUE)
write.csv(boundMatrixFix, paste(projectPath, "Rfiles/boundMatrix.csv", sep = ""), row.names = FALSE)
