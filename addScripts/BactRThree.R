library(dplyr)
library(ape)
library(reshape2)
library(stringr)
library(pheatmap)
library(phangorn)
library(dendextend)
library(argparse)

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
            write.csv(repTrfIds, paste(resultDIR, "repTrfIds.csv", sep=""), row.names=TRUE)
          }else{write("repTRFIds lacking data, unable to create file.", file = err, append = TRUE)}
        }else{write("Lacking data in trf txt, unable to create repTRFIds.csv.", file = err, append = TRUE)}
      }else{write("Missing trfParsed.txt, unable to create repTRFIds.csv.", file = err, append = TRUE)}
    }else{write("Lacking data in repeat csv, unable to create repTRFIds.csv.", file = err, append = TRUE)}
  }else{write("Missing RepeatNames.csv, unable to create repTRFIds.csv.", file = err, append = TRUE)}
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
        write.FASTA(rvdFASTA,paste(resultDIR, "rvdIds.FASTA", sep=""))
      }else{write("rvdFASTA in rvd_ids_file function lacking data, unable to create file", file = err, append = TRUE)}
    }else{write("rvdFASTA lacking data, unable to create file", file = err, append = TRUE)}
  }else{write("rvdFASTA was not generated, unable to modify file", file = err, append = TRUE)}
}

rvd_nucs_file <- function()
{
  if(file.exists(rvdNucsFile))
  {
    rvdNuc <- read.csv(rvdNucsFile)
    if(object.size(rvdNucs) > 0)
    {
      newRvd <- rvdNucs[rvdNucs$ID %in% talIDS,]
      write.csv(newRvd,paste(resultDIR, "rvdIDs.csv"))
    }else{write("rvdNucs lacking data, unable to modify file", file = err, append = TRUE)}
  }else{write("rvdNucs.csv was not generated, unable to modify file", file = err, append = TRUE)}
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
            write.FASTA(faaOrthoOrder,paste(resultDIR, "faaIds.FASTA", sep=""))
          }else{write("faaOrthoOrder in ortho_melt function lacking data, unable to create file", file = err, append = TRUE)}
        }else{write("faaFASTA file lacking data, unable to continue with faaIds.FASTA", file = err, append = TRUE)}
      }else{write("Missing faaFASTA, unable to continue with faaIds.FASTA", file = err, append = TRUE)}
    }else{write("Orthogroups txt lacking data, unable to continue with faaIds.FASTA", file = err, append = TRUE)}
  }else{write("Missing Orthogroups.txt, unable to continue with faaIds.FASTA", file = err, append = TRUE)}
}

#Open error file for writing
errFile <- paste(Logging, "RErrors.txt", sep = "")
err <- file(errFile)

#Creating ideas
if(file.exists(scoaryCSVFile))
{
  scoaryCSV <- read.csv(scoaryCSVFile)
  if(nrow(scoaryCSV) > 0)
  {
    ids <- scoaryCSV[scoaryCSV$Bonferroni_p<0.05,1]
    rep_trf_file()
    ortho_melt()
    if(file.exists(talGroupsCSVFile))
    {
      talGroupsCSV <- read.csv(talGroupsCSVFile, stringsAsFactors = FALSE)
      if(nrow(talGroupsCSV) > 0)
      {
        talIDS <- talGroupsCSV$TAL[paste("TALGroup_",talGroupsCSV$Group,sep="") %in% ids]
        rvd_ids_file()
        rvd_nucs_file()
      }else{write("talGroups csv lacking data, unable to create talIDS", file = err, append = TRUE)}
    }else{write("Missing talgroups csv, unable to create talIDS", file = err, append = TRUE)}
  }else{write("Scoary csv lacking data, unable to create ids", file = err, append = TRUE)}
}else{write("Missing scoary csv, unable to create ids", file = err, append = TRUE)}
close(err)
