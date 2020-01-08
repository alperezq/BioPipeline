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
orthoGCMatrixFile <- paste(projectPath, "Rfiles/OrthoGCMatrix.csv", sep="")
matTRFFile <- paste(projectPath, "Rfiles/matTRF.csv", sep="")
groupMatFile <- paste(projectPath, "Rfiles/distal_GroupMatrix.csv", sep="")

#Checks existance of files, converts to data frames and binds them, manipulates, writes Bound Matrix
boundMatrixCreation <- function()
{
  #Counter to track items in list, list to hold working dataframes
  i <- 0
  dataList <- list()

  if(file.exists(orthoGCMatrixFile)){
    orthoGCMatrix <- read.csv(orthoGCMatrixFile)
    if(object.size(orthoGCMatrix) > 0){
      dataList[[i]] <- orthoGCMatrix
      i <- i + 1
    }else{print("OrthoGCMatrix.csv was generated, but has no data.")}
  }else{print("OrthoGCMatrix.csv was not generated, not using for Bound Matrix.")}

  if(file.exists(matTRFFile)){
    matTRF <- read.csv(matTRFFile)
    if(object.size(matTRF) > 0){
      dataList[[i]] <- matTRF
      i <- i + 1
    }else{print("matTRF.csv was generated, but has no data.")}
  }else{print("matTRF.csv was not generated, not using for Bound Matrix.")}

  if(file.exists(groupMatFile)){
    groupMat <- read.csv(groupMatFile)
    if(object.size(groupMat) > 0){
      dataList[[i]] <- groupMatFile
      i <- i + 1
    }else{print("distal_GroupMatrix.csv was generated, but has no data.")}
  }else{print("distal_GroupMatrix.csv was not generated, not using for Bound Matrix.")}

  #Creation of bound matrix after manipulations in python3
  if(length(dataList) > 0){
    boundMatrix <- bind_rows(dataList)
    boundMatrixFix <- slice(boundMatrix, 1:(n()-2))
    colnames(boundMatrixFix) <- gsub(".","-", colnames(boundMatrixFix), fixed = TRUE)
    write.csv(boundMatrixFix, paste(projectPath, "Rfiles/boundMatrix.csv", sep = ""), row.names = FALSE)
  }else{print("None of required csv's are available, unable to generate boundMatrix.csv.")}
}

#Call to above function
boundMatrixCreation()
