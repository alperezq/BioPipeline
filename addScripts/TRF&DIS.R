#!/usr/bin/env RScript
#R script for first section of CSU Bioinformatics Pipeline

#Use Pacman for checking and loading of other libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2,pheatmap, plyr, ComplexHeatmap, reshape2, ape, dplyr, stringr, argparse)

#Use of argparse to create parser to read in project path
args <- commandArgs(trailingOnly = TRUE)

projectPath = args[1]

TRFPath = paste(projectPath, "TRFfiles/", sep = "")
DisTALPath = paste(projectPath, "DISTALfiles/", sep = "")
RPath = paste(projectPath, "Rfiles/", sep = "")


repeatsTRF <- read.delim(paste(TRFPath,"trfParsed.txt", sep = ""), header = FALSE, sep = " ")

disOut <- as.matrix(read.table(paste(DisTALPath, "disTALOut.mat", sep = ""), header=TRUE, sep = "\t",row.names = 1, as.is=TRUE))


countsTRF<-as.data.frame(table(repeatsTRF$V6,repeatsTRF$V1))
matTRF<-dcast(countsTRF, Var1 ~ Var2, value.var = "Freq")
aggTRF<-repeatsTRF[,c(1,5,6)]
repSumsTRF<-as.data.frame(aggregate(aggTRF$V5 ~ aggTRF$V6 + aggTRF$V1, FUN= sum))
colnames(repSumsTRF)<-c("V1","V2","V3")
matTRF<-dcast(repSumsTRF, V1 ~ V2, value.var = "V3")
matTRF[is.na(matTRF)]<-0

pdf(paste(RPath, "matTRFPheatmap.pdf", sep = ""))
pheatmap(matTRF[,2:ncol(matTRF)])
dev.off()


tr <- nj((disOut))
bla<-disOut[,1:(ncol(disOut))-1]
a <-hclust(as.dist(bla)) #Cut tree and get groups
CUT<-as.data.frame(cutree(a,h=0.5))
Groups <- as.data.frame(cbind(rownames(disOut),CUT[,1]))
colnames(Groups)<-c("TAL","Group")
write.csv(Groups, paste(RPath, "disOutR.csv", sep = ""), row.names = TRUE)
