#!/usr/bin/env RScript
#R script for first section of CSU Bioinformatics Pipeline

#Use Pacman for checking and loading of other libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2,pheatmap, plyr, ComplexHeatmap, reshape2, ape, dplyr, stringr, argparse)

#Use of argparse to create parser to read in project path
args <- commandArgs(trailingOnly = TRUE)

projectPath = args[1]
projectPath = "/Users/tyr/Documents/CSU/Plant/TEST/"

orthoPath = paste(projectPath, "ORTHOfiles/", sep = "")
TRFPath = paste(projectPath, "TRFfiles/", sep = "")
DisTALPath = paste(projectPath, "DISTALfiles/", sep = "")
RPath = paste(projectPath, "Rfiles/", sep = "")

orthoGC <- read.delim(paste(orthoPath,"Orthogroups.GeneCount.csv", sep=""))

repeatsTRF <- read.delim(paste(TRFPath,"trfParsed.txt", sep = ""), header = FALSE, sep = " ")

disOut <- as.matrix(read.table(paste(DisTALPath, "disTALOut.mat", sep = ""), header=TRUE, sep = "\t",row.names = 1, as.is=TRUE))


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

Groups <- read.delim(paste(DisTALPath, "disTALOut.TALgroups.txt", sep =""))
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

write.csv(boundMatrixFix, paste(RPath, "boundMatrix.csv", sep = ""), row.names = TRUE)


