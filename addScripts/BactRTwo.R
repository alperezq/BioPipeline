#!/usr/bin/env Rscript
#Declaration of libraries
library(dplyr)
library(argparse)

#Take in args from IdunsBridge
args <- commandArgs(trailingOnly=TRUE)
projectPath = args[1]

#Assign paths to variables
orthoGCMatrixFile <- paste(projectPath, "Rfiles/OrthoGCMatrix.csv", sep="")
matTRFFile <- paste(projectPath, "Rfiles/matTRF.csv", sep="")
groupMatFile <- paste(projectPath, "Rfiles/distal_GroupMatrix.csv", sep="")

#Checks existance of files, converts to data frames and binds them, manipulates, writes Bound Matrix
boundMatrixCreation <- function()
{
  #Counter to track items in list, list to hold working dataframes
  i <- 0

  orthoGCMatrix <- data.frame()
  matTRF <- data.frame()
  groupMat <- data.frame()

  if(file.exists(orthoGCMatrixFile)){
    orthoGCMatrix <- read.csv(orthoGCMatrixFile)
    if(object.size(orthoGCMatrix) > 0){
      i <- i + 1
    }else{write("OrthoGCMatrix.csv was generated, but has no data.", file = err, append = TRUE)}
  }else{write("OrthoGCMatrix.csv was not generated, not using for Bound Matrix.", file = err, append = TRUE)}

  if(file.exists(matTRFFile)){
    matTRF <- read.csv(matTRFFile)
    if(object.size(matTRF) > 0){
      i <- i + 1
    }else{write("matTRF.csv was generated, but has no data.", file = err, append = TRUE)}
  }else{write("matTRF.csv was not generated, not using for Bound Matrix.", file = err, append = TRUE)}

  if(file.exists(groupMatFile)){
    groupMat <- read.csv(groupMatFile)
    if(object.size(groupMat) > 0){
      i <- i + 1
    }else{write("distal_GroupMatrix.csv was generated, but has no data.", file = err, append = TRUE)}
  }else{write("distal_GroupMatrix.csv was not generated, not using for Bound Matrix.", file = err, append = TRUE)}

  #Creation of bound matrix after manipulations in python3
  if(i > 0){
    boundMatrix <- bind_rows(orthoGCMatrix, matTRF, groupMat)
    boundMatrixFix <- slice(boundMatrix, 1:(n()-2))
    colnames(boundMatrixFix) <- gsub(".","-", colnames(boundMatrixFix), fixed = TRUE)
    write.csv(boundMatrixFix, paste(projectPath, "Rfiles/boundMatrix.csv", sep = ""), row.names = FALSE)
  }else{write("None of required csv's are available, unable to generate boundMatrix.csv.", file = err, append = TRUE)}
}

#Call to above function
errFile <- paste(projectPath, "Rfiles/RErrors.txt", sep = "")
err <- file(errFile)
boundMatrixCreation()
close(err)
