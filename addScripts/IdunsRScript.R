#!/usr/bin/env RScript
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

#Establish paths based on project presets,
orthoPath = paste(projectPath, "ORTHOfiles/", sep = "")
RPath = paste(projectPath, "Rfiles/", sep = "")
TRFPath = paste(projectPath, "TRFfiles/", sep = "")
DisTALPath = paste(projectPath, "DISTALfiles/", sep = "")

#Assign data files to variables
orthoGC <- read.delim(paste(orthoPath,"Orthogroups.GeneCount.csv", sep=""))
repeatsTRF <- read.delim(paste(TRFPath,"trfParsed.txt", sep = ""), header = FALSE, sep = " ")
disOut <- as.matrix(read.table(paste(DisTALPath, "Outputs/disTALOut.mat", sep = ""), header=TRUE, sep = "\t",row.names = 1, as.is=TRUE))
statsPerSpecies <- read.delim(paste(orthoPath,"Statistics_PerSpecies.csv", sep=""))
Groups <- read.delim(paste(DisTALPath, "Outputs/disTALOut.TALgroups.csv", sep =""))


#Manipulation of orthoGC to Matrix
orthoGCMatrix <- as.matrix(orthoGC[,2:(ncol(orthoGC)-1)], rownames.force = 0, nrow(5), ncol(5))
orthoGCMatrix[orthoGCMatrix>1]<-1
Z1Orth <- orthoGC[,1:(ncol(orthoGC)-1)]
Z1Orth[,2:ncol(Z1Orth)] <-apply(Z1Orth[,2:ncol(Z1Orth)], 2, function(x) ifelse(x > 1, 1, x))

#statsPerSpecies manipulations and creation
statsPerc <- statsPerSpecies[statsPerSpecies$X == "Percentage of genes in species-specific orthogroups",]
statsPerc <- as.data.frame(t(statsPerc))
write.csv(statsPerc, paste(RPath, "statsPerSpeciesR.csv", sep = ""),  row.names = TRUE)

#Creation of various plots using OrthoGC data
pdf(paste(RPath, "orthoGCggplot.pdf", sep = ""))
ggplot(data=orthoGC,aes(x=AUST2013, y=JW11089)) + geom_point(size=1, shape=6)
dev.off()
pdf(paste(RPath, "orthoGCgg&geom.pdf", sep = ""))
ggplot(data=orthoGC,aes(x=AUST2013, y=JW11089)) + geom_point()
dev.off()
pdf(paste(RPath,"orthoGCMatrixHeatmap.pdf", sep = ""))
Heatmap(orthoGCMatrix, row_order = sort(rownames(orthoGCMatrix)),)
dev.off()
pdf(paste(RPath,"orthoGCPheatmap.pdf", sep = ""))
pheatmap(orthoGC[,2:(ncol(orthoGC)-1)])
dev.off()
pdf(paste(RPath, "orthoGCPheatmap.pdf", sep = ""))
pheatmap(orthoGC[,2:(ncol(orthoGC)-1)])
dev.off()
pdf(paste(RPath, "orthoGCBoxplot.pdf", sep = ""))
boxplot(orthoGC[,2:(ncol(orthoGC)-1)])
dev.off()
pdf(paste(RPath, "Z10rthPheatmap.pdf", sep = ""))
pheatmap(Z1Orth[,2:ncol(Z1Orth)])
dev.off()


#TandemRepeatFinder data manipulations
aggTRF<-repeatsTRF[,c(1,5,7)]
repSumsTRF<-as.data.frame(aggregate(aggTRF$V5 ~ aggTRF$V7 + aggTRF$V1, FUN= sum))
colnames(repSumsTRF)<-c("V1","V2","V3")
matTRF<-dcast(repSumsTRF, V1 ~ V2, value.var = "V3")
matTRF[is.na(matTRF)]<-0

RepNames<-data.frame("Sequence"=matTRF$V1)
RepNames$Name<-paste("Repeats",1:nrow(RepNames),sep="_")
matTRF$V1<-RepNames$Name
write.csv(matTRF, paste(RPath, "matTRF.csv", sep = ""), row.names = TRUE)
write.csv(RepNames, paste(RPath, "RepeatNames.csv", sep = ""), row.names = TRUE)

pdf(paste(RPath, "matTRFPheatmap.pdf", sep = ""))
pheatmap(matTRF[,2:ncol(matTRF)])
dev.off()


#DisTAL data manipulations
Groups$Group<-gsub("^.* ","",Groups$TAL.Group)
Groups$Genome<-str_split_fixed(Groups$TAL.Group,"\\|",2)[,1]
Groups$Group<-paste("TALGroup_",Groups$Group,sep="")
G<-Groups[,2:3]
GroupMat<-dcast(G,Group ~ Genome)
write.csv(GroupMat, paste(RPath, "distal_GroupMatrix.csv", sep = ""), row.names = TRUE)


#Formatting and gathering of matrixes for creation of bound matrix
orthoGCMatrix <- as.matrix(orthoGC[,1:(ncol(orthoGC)-1)], rownames.force = 0, nrow(5), ncol(5))
colnames(orthoGCMatrix) <- gsub(pattern = "\\.", "-", x = colnames(orthoGCMatrix))
colnames(GroupMat)[colnames(GroupMat)=="Group"] <- "V1"
colnames(orthoGCMatrix)[colnames(orthoGCMatrix)=="X"] <- "V1"
boundMatrix <- rbind(orthoGCMatrix, matTRF, GroupMat, sep = ",", stringsAsFactors = FALSE)
boundMatrixFix <- slice(boundMatrix, 1:(n()-1))
write.csv(boundMatrixFix, paste(RPath, "boundMatrix.csv", sep = ""), row.names = FALSE)
