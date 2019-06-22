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
RPath = paste(projectPath, "Rfiles/", sep = "")


orthoGC <- read.delim(paste(orthoPath,"Orthogroups.GeneCount.csv", sep=""))

statsPerSpecies <- read.delim(paste(orthoPath,"Statistics_PerSpecies.csv", sep=""))


orthoGCMatrix <- as.matrix(orthoGC[,2:(ncol(orthoGC)-1)], rownames.force = 0, nrow(5), ncol(5))
orthoGCMatrix[orthoGCMatrix>1]<-1

Z1Orth <- orthoGC[,1:(ncol(orthoGC)-1)]
Z1Orth[,2:ncol(Z1Orth)] <-apply(Z1Orth[,2:ncol(Z1Orth)], 2, function(x) ifelse(x > 1, 1, x))


statsPerc <- statsPerSpecies[statsPerSpecies$X == "Percentage of genes in species-specific orthogroups",]
statsPerc <- as.data.frame(t(statsPerc))
write.csv(statsPerc, paste(RPath, "statsPerSpeciesR.csv", sep = ""),  row.names = TRUE)


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
