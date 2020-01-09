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

scoaryCSV = read.csv(args[1])
repeatCSV = read.csv(args[2])
boundCSV = read.csv(args[3])
traitCSV = read.csv(args[4])
talGroupsCSV = read.csv(args[5], stringsAsFactors = FALSE)
trfTXT = read.delim(args[6],sep = " ",header = FALSE)
orthogroupsTXT = read.delim(args[7],sep=" ",header = FALSE)
kTree = read.tree(args[8])
rvdFASTA = read.FASTA(args[9],type="AA")
resultDIR = args[10]
conCSV = read.csv(args[11])
faaFASTA = read.FASTA(args[12],type="AA")
rvdNucs = read.csv(args[13])


#Functions for creation of datasets
rep_trf_file <- function()
{
  if(object.size(repeatCSV) > 0)
  {
    if(object.size(trfTXT) > 0)
    {
      repeatWithTrf <- left_join(repeatCSV,trfTXT,by=c("Sequence"="V7"))
      repTrfIds <- repeatWithTrf[repeatWithTrf$Name %in% ids,]
      if(nrow(repTrfIds) > 0)
      {
        write.csv(repTrfIds, paste(resultDIR, "repTrfIds.csv", sep=""), row.names=TRUE)
      }else{write("repTRFIds lacking data, unable to create file", file = err, append = TRUE)}
    }else{write("Lacking data in trf txt, unable to create repTRFIds csv", file = err, append = TRUE)}
  }else{write("Lacking data in repeat csv, unable to create repTRFIds csv", file = err, append = TRUE)}
}

rvd_ids_file <- function()
{
  if(object.size(rvdFASTA) > 0)
  {
    rvdFASTA <- rvdFASTA[names(rvdFASTA) %in% talIDS]
    if(object.size(rvdFASTA) > 0)
    {
      write.FASTA(rvdFASTA,paste(resultDIR, "rvdIds.FASTA", sep=""))
    }else{write("rvdFASTA in rvd_ids_file function lacking data, unable to create file", file = err, append = TRUE)}
  }else{write("rvdFASTA lacking data, unable to create file", file = err, append = TRUE)}
}

rvd_nucs_file <- function()
{
  if(object.size(rvdNucs) > 0)
  {
    newRvd <- rvdNucs[rvdNucs$ID %in% talIDS,]
    write.csv(newRvd,paste(resultDIR, "rvdIDs.csv"))
  }else{write("rvdNucs lacking data, unable to modify file", file = err, append = TRUE)}
}

ortho_melt <- function()
{
  if(object.size(orthogroupsTXT) > 0)
  {
    orthoMelt <- melt(orthogroupsTXT,id.vars = "V1")
    orthoMelt <- orthoMelt[,c(1,3)]
    orthoMelt <- orthoMelt[orthoMelt$value != "",]
    orthoMelt$V1 <- gsub(":", "",orthoMelt$V1)

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
    }else{write("faaFASTA file lacking data, unable to continue with Ortho Melt function", file = err, append = TRUE)}
  }else{write("Orthogroups txt lacking data, unable to continue with Ortho Melt function", file = err, append = TRUE)}
}

#Open error file for writing
errFile <- paste(projectPath, "Rfiles/RErrors.txt", sep = "")
err <- file(errFile)

#Creating ideas
if(nrow(scoaryCSV) > 0)
{
  ids <- scoaryCSV[scoaryCSV$Bonferroni_p<0.05,1]
  rep_trf_file()
  ortho_melt()
  rvd_nucs_file()
  if(nrow(talGroupsCSV) > 0)
  {
    talIDS <- talGroupsCSV$TAL[paste("TALGroup_",talGroupsCSV$Group,sep="") %in% ids]
    rvd_ids_file()
  }else{write("talGroups csv lacking data, unable to create talIDS", file = err, append = TRUE)}
}else{write("Scoary csv lacking data, unable to create ids", file = err, append = TRUE)}

close(err)
