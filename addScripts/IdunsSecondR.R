library(dplyr)
library(ape)
library(reshape2)
library(stringr)
library(pheatmap)
library(phangorn)
library(dendextend)
library(argparse)
#541-828896

#collects arguments
args <- commandArgs(trailingOnly = TRUE)

scoaryCSV = read.csv("Downloads/idunsTEST/SCOARYfiles/Race_9b_14_07_2019_1415.results.csv")
repeatCSV = read.csv("Downloads/idunsTEST/Rfiles/RepeatNames.csv")
boundCSV = read.csv("Downloads/idunsTEST/Rfiles/boundMatrix.csv")
traitCSV = read.csv("Downloads/Trait.csv")
talGroupsCSV = read.csv("Downloads/idunsTEST/DISTALfiles/Outputs/disTALOut.TALgroups.csv", stringsAsFactors = FALSE)
trfTXT = read.delim("Downloads/idunsTEST/TRFfiles/trfParsed.txt",sep = " ",header = FALSE)
orthogroupsTXT = read.delim("Downloads/idunsTEST/ORTHOfiles/Orthogroups.txt",sep=" ",header = FALSE)
kTree = read.tree("Downloads/idunsTEST/KSNP3files/kSNP3_results/tree.ML.tre")
rvdFASTA = read.FASTA("Downloads/idunsTEST/DISTALfiles/rvdCombo.FASTA")
resultDIR = "Downloads/idunsTEST/Results/"
conCSV = read.csv("Downloads/Concatenated_createdbyalvaro.csv")
faaFASTA = read.FASTA("Downloads/AllFAAs_withnames.fa")

#Gathering of file/dir paths
scoaryCSV <- read.csv(args[1])
repeatCSV <- read.csv(args[2])
boundCSV <- read.csv(args[3])
traitCSV <- read.csv(args[4])
talGroupsCSV <- read.csv(args[5])
trfTXT <- read.delim(args[6],sep = " ",header = FALSE)
orthogroupsTXT <- read.delim(args[7],sep=" ",header = FALSE)
kTree <- read.tree(args[8])
rvdFASTA <- read.FASTA(args[9])
resultDIR = args[10]
conCSV <- read.delim(args[11],sep="\t")
faaFASTA <- read.FASTA(args[12])

#Creating Ids
ids <- scoaryCSV[scoaryCSV$Bonferroni_p<0.05,1]
talIDS <- talGroupsCSV$TAL[paste("TALGroup_",talGroupsCSV$Group,sep="") %in% ids]

#Manipulation of repeats file
repeatWithTrf <- left_join(repeatCSV,trfTXT,by=c("Sequence"="V7"))
repTrfIds <- repeatWithTrf[repeatWithTrf$Name %in% ids,]
write.csv(repTrfIds, paste(resultDIR, "repTrfIds.csv", sep=""), row.names=TRUE)

#Manipulation of concatenated file
rvdFASTA <- rvdFASTA[names(rvdFASTA) %in% talIDS]
write.FASTA(rvdFASTA,paste(resultDIR, "rvdIds.FASTA", sep=""))

#Manipulation of Orthogroups.txt
orthoMelt <- melt(orthogroupdsTXT,id.vars = "V1")
orthoMelt <- orthoMelt[,c(1,3)]
orthoMelt <- orthoMelt[orthoMelt$value != "",]
orthoMelt$V1 <- gsub(":", "",orthoMelt$V1)

#Manipulation of FAA FASTA, combine with Ortho
fastaIds <-data.frame("Short_ID"=str_split_fixed(names(faaFASTA)," ",3)[,2],"fastaID"=names(faaFASTA))
orthoFastaIds <- left_join(fastaIds,orthoMelt,by=c("Short_ID"="value"))
orthoFastaIds$New_ID<-paste(orthoFastaIds$V1,orthoFastaIds$fastaID,sep=":")
names(faaFASTA) <- orthoFastaIds$New_ID
orthoSort <- sort(orthoFastaIds$New_ID[orthoFastaIds$V1 %in% ids])
faaOrthoSort <- faaFASTA[names(faaFASTA) %in% orthoSort]
faaOrthoOrder <- faaOrthoSort[order(names(faaOrthoSort))]
write.FASTA(faaOrthoOrder,paste(resultDIR, "faaIds.FASTA", sep=""))
