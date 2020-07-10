#!/usr/bin/env Rscript
#Declaration of libraries
library(reshape2)
library(dplyr)
library(argparse)
library(ape)
library(stringr)

#collects arguments
args <- commandArgs(trailingOnly = TRUE)

scoaryCSVFile <- args[1]
repeatCSVFile <- args[2]


talGroupsCSVFile <- args[3]
trfTXTFile <- args[4]
orthogroupsTXTFile <- args[5]

rvdFASTAFile <- args[6]
resultDIR <- args[7]

faaFASTAFile <- args[8]
rvdNucsFile <- args[9]
Logging <- args[10]

#Functions for creation of datasets
rep_trf_file <- function()
{

  if(file.exists(repeatCSVFile))
  {
    repeatCSV <- read.csv(repeatCSVFile)
    if(object.size(repeatCSV) > 0)
    {
      if(file.exists(trfTXTFile))
      {
        trfTXT <- read.delim(trfTXTFile, sep = " ", header = FALSE)
        if(object.size(trfTXT) > 0)
        {
          repeatWithTrf <- left_join(repeatCSV,trfTXT,by=c("Sequence"="V7"))
          repTrfIds <- repeatWithTrf[repeatWithTrf$Name %in% ids,]
          if(nrow(repTrfIds) > 0)
          {
            write.csv(repTrfIds, paste(resultDIR, "sc_repTrfIds.csv", sep=""), row.names=TRUE)
          }else{write("BactRThree: repTRFIds lacking data, unable to create file.", file = err, append = TRUE)}
        }else{write("BactRThree: Lacking data in trf txt, unable to create repTRFIds.csv.", file = err, append = TRUE)}
      }else{write("BactRThree: Missing trfParsed.txt, unable to create repTRFIds.csv.", file = err, append = TRUE)}
    }else{write("BactRThree: Lacking data in repeat csv, unable to create repTRFIds.csv.", file = err, append = TRUE)}
  }else{write("BactRThree: Missing RepeatNames.csv, unable to create repTRFIds.csv.", file = err, append = TRUE)}
}

ortho_melt <- function()
{
  if(file.exists(orthogroupsTXTFile))
  {
    orthogroupsTXT <- read.delim(orthogroupsTXTFile,sep=" ",header = FALSE)
    if(object.size(orthogroupsTXT) > 0)
    {
      orthoMelt <- melt(orthogroupsTXT,id.vars = "V1")
      orthoMelt <- orthoMelt[,c(1,3)]
      orthoMelt <- orthoMelt[orthoMelt$value != "",]
      orthoMelt$V1 <- gsub(":", "",orthoMelt$V1)

      if(file.exists(faaFASTAFile))
      {
        faaFASTA <- read.FASTA(faaFASTAFile,type="AA")
        if(object.size(faaFASTA) > 0)
        {
          fastaIds <-data.frame("Short_ID"=str_split_fixed(names(faaFASTA)," ",3)[,2],"fastaID"=names(faaFASTA))
          orthoFastaIds <- left_join(fastaIds,orthoMelt,by=c("Short_ID"="value"))
          orthoFastaIds$New_ID <- paste(orthoFastaIds$V1,orthoFastaIds$fastaID,sep=":")
          names(faaFASTA) <- orthoFastaIds$New_ID
          orthoSort <- sort(orthoFastaIds$New_ID[orthoFastaIds$V1 %in% ids])
          faaOrthoSort <- faaFASTA[names(faaFASTA) %in% orthoSort]
          faaOrthoOrder <- faaOrthoSort[order(names(faaOrthoSort))]

          if(object.size(faaOrthoOrder) > 0)
          {
            write.FASTA(faaOrthoOrder,paste(resultDIR, "sc_faaIds.FASTA", sep=""))
          }else{write("BactRThree: faaOrthoOrder in ortho_melt function lacking data, unable to create file", file = err, append = TRUE)}
        }else{write("BactRThree: faaFASTA file lacking data, unable to continue with faaIds.FASTA", file = err, append = TRUE)}
      }else{write("BactRThree: Missing faaFASTA, unable to continue with faaIds.FASTA", file = err, append = TRUE)}
    }else{write("BactRThree: Orthogroups txt lacking data, unable to continue with faaIds.FASTA", file = err, append = TRUE)}
  }else{write("BactRThree: Missing Orthogroups.txt, unable to continue with faaIds.FASTA", file = err, append = TRUE)}
}

rvd_ids_file <- function()
{
  if(file.exists(rvdFASTAFile))
  {
    rvdFASTA <- read.FASTA(rvdFASTAFile,type="AA")
    if(object.size(rvdFASTA) > 0)
    {
      rvdFASTA <- rvdFASTA[names(rvdFASTA) %in% talIDS]
      if(object.size(rvdFASTA) > 0)
      {
        write.FASTA(rvdFASTA,paste(resultDIR, "sc_rvdIDs.FASTA", sep=""))
      }else{write("BactRThree: rvdFASTA in rvd_ids_file function lacking data, unable to create file", file = err, append = TRUE)}
    }else{write("BactRThree: rvdFASTA lacking data, unable to create file", file = err, append = TRUE)}
  }else{write("BactRThree: rvdFASTA was not generated, unable to modify file", file = err, append = TRUE)}
}

rvd_nucs_file <- function()
{
  if(file.exists(rvdNucsFile))
  {
    rvdNucs <- read.csv(rvdNucsFile)
    if(object.size(rvdNucs) > 0)
    {
      newRvd <- rvdNucs[rvdNucs$ID %in% talIDS,]
      write.csv(newRvd,paste(resultDIR, "sc_rvdIDs.csv", sep=""))
    }else{write("BactRThree: rvdNucs lacking data, unable to modify file", file = err, append = TRUE)}
  }else{write("BactRThree: rvdNucs.csv was not generated, unable to modify file", file = err, append = TRUE)}
}

#Open error file for writing
errFile <- paste(Logging, "RErrorsBact3.txt", sep = "")
err <- file(errFile)

#Creating ids
if(file.exists(scoaryCSVFile))
{
	#Change to read delim for bayes results
  scoaryCSV <- read.csv(scoaryCSVFile)
  if(nrow(scoaryCSV) > 0)
  {
	   #Filter out values based on Bonferroni_p
    ids <- as.vector(scoaryCSV[scoaryCSV$Bonferroni_p<0.05,1])
    if(length(ids) > 0)
    {
      rep_trf_file()
      ortho_melt()
      if(file.exists(talGroupsCSVFile))
      {
        talGroupsCSV <- read.csv(talGroupsCSVFile, stringsAsFactors = FALSE)
        if(nrow(talGroupsCSV) > 0)
        {
          talIDS <- as.vector(talGroupsCSV$TAL[paste("TALGroup_",talGroupsCSV$Group,sep="") %in% ids])
          if(length(talIDS) > 0){
            rvd_ids_file()
            rvd_nucs_file()
          }else{write("BactRThree: talIDS has no data, unable to create talIDS", file = err, append = TRUE)}
        }else{write("BactRThree: talGroups csv lacking data, unable to create talIDS", file = err, append = TRUE)}
      }else{write("BactRThree: Missing talgroups csv, unable to create talIDS", file = err, append = TRUE)}
    }else{write("BactRThree: ids has no data, unable to continue processing heres", file = err, append = TRUE)}
  }else{write("BactRThree: Scoary csv lacking data, unable to create ids", file = err, append = TRUE)}
}else{write("BactRThree: Missing scoary csv, unable to create ids", file = err, append = TRUE)}
close(err)
