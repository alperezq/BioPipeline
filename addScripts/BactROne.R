#!/usr/bin/env Rscript
#R script for first section of CSU Bioinformatics Pipeline

#Declaration of libraries
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(reshape2)
library(ape)
library(stringr)
library(plyr)
library(dplyr)
library(argparse)

#Use of argparse to create parser to read in project path
args <- commandArgs(trailingOnly = TRUE)
projectPath = args[1]

#Assign data files to variables
orthoGCFile <- paste(projectPath,"ORTHOfiles/Orthogroups.GeneCount.csv", sep="")
repeatsTRFFile <- paste(projectPath, "TRFfiles/trfParsed.txt", sep = "")
statsPerSpeciesFile <- paste(projectPath,"ORTHOfiles/Statistics_PerSpecies.csv", sep = "")
GroupsFile <- paste(projectPath, "DISTALfiles/Outputs/disTALOut.TALgroups.csv", sep ="")

#Manipulation of orthoGC to Matrix
OrthoGC <- function()
{
  if(file.exists(orthoGCFile)){
    orthoGC <- read.delim(orthoGCFile)
    if(object.size(orthoGC) > 0)
    {
      orthoGCMatrix <- as.matrix(orthoGC[,1:(ncol(orthoGC)-1)], rownames.force = 0, nrow(5), ncol(5))
      colnames(orthoGCMatrix) <- gsub(pattern = "\\.", "-", x = colnames(orthoGCMatrix))
      colnames(orthoGCMatrix)[colnames(orthoGCMatrix)=="X"] <- "V1"
      write.csv(orthoGCMatrix, paste(RPath, "OrthoGCMatrix.csv", sep=""), row.names = FALSE)
    }else{write("Orthogroups.GeneCount.csv was generated, but has no data", file = err, append = TRUE)}
  }else{write("Orthogroups.GeneCount.csv was not generated.", file = err, append = TRUE)}
}

#Manipulation of statsPerSpecies
perSpeciesStats <- function()
{
  if(file.exists(statsPerSpeciesFile)){
    statsPerSpecies <- read.delim(statsPerSpeciesFile)
    if(object.size(statsPerSpecies) > 0){
      statsPerc <- statsPerSpecies[statsPerSpecies$X == "Percentage of genes in species-specific orthogroups",]
      statsPerc <- as.data.frame(t(statsPerc))
      write.csv(statsPerc, paste(RPath, "statsPerSpeciesR.csv", sep = ""),  row.names = TRUE)
    }else{write("Statistics_PerSpecies.csv was generated, but has no data.", file = err, append = TRUE)}
  }else{write("Statistics_PerSpecies.csv was not generated", file = err, append = TRUE)}
}

#Manipulation of trf -> matTRf and repeatNames
TandemRepeatFinder <- function()
{
  if(file.exists(repeatsTRFFile)){
    repeatsTRF <- read.delim(repeatsTRFFile, header = FALSE, sep = "")
    if(object.size(repeatsTRF) > 0){
      aggTRF<-repeatsTRF[,c(1,5,7)]
      repSumsTRF<-as.data.frame(aggregate(aggTRF$V5 ~ aggTRF$V7 + aggTRF$V1, FUN= sum))
      colnames(repSumsTRF)<-c("V1","V2","V3")
      matTRF<-dcast(repSumsTRF, V1 ~ V2, value.var = "V3")
      matTRF[is.na(matTRF)]<-0

      RepNames<-data.frame("Sequence"=matTRF$V1)
      RepNames$Name<-paste("Repeats",1:nrow(RepNames),sep="_")
      matTRF$V1<-RepNames$Name
      colnames(matTRF) <- gsub(pattern = "\\.", "-", x = colnames(matTRF))
      write.csv(matTRF, paste(RPath, "matTRF.csv", sep = ""), row.names = FALSE)
      write.csv(RepNames, paste(RPath, "RepeatNames.csv", sep = ""), row.names = FALSE)
    }else{write("trfParsed.txt was generated, but has no data.", file = err, append = TRUE)}
  }else{write("trfParsed.txt as not generated", file = err, append = TRUE)}
}

#manipulation of Groups to Talgroups.csv
DisTAL <- function()
{
  if(file.exists(GroupsFile)){
    Groups <- read.csv(GroupsFile)
    if(object.size(Groups) > 0){
      Groups$Genome<-str_split_fixed(Groups$TAL,"\\|",2)[,1]
      Groups$Group<-paste("TALGroup_",Groups$Group,sep="")
      G<-Groups[,2:3]
      GroupMat<-dcast(G,Group ~ Genome)
      colnames(GroupMat) <- gsub(pattern = "\\.", "-", x = colnames(GroupMat))
      colnames(GroupMat)[colnames(GroupMat)=="Group"] <- "V1"
      write.csv(GroupMat, paste(RPath, "distal_GroupMatrix.csv", sep = ""), row.names = FALSE)
    }else{write("disTALOut.TALgroups.csv was generated, but has no data", file = err, append = TRUE)}
  }else{write("disTALOut.TALgroups.csv was not generated.", file = err, append = TRUE)}
}


#Creation of text file for errors
errFile <- paste(projectPath, "Rfiles/RErrors.txt", sep = "")
file.create(errFile)
err <- file(errFile)

#Calls to above functions
OrthoGC()
perSpeciesStats()
TandemRepeatFinder()
DisTAL()

close(err)
