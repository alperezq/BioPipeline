library(ape)
library(pheatmap)
library(stringr)
library(dplyr)
library(phangorn)
library(argparse)

args <- commandArgs(trailingOnly=TRUE)
projectPath <- args[1]

######## read and modify phylogenetic tree
Ktree <- read.tree(paste(projectPath, "KSNP3files/kSNP3_results/tree.ML.tre", sep=""))

Ktree$node.label<-NULL #delete node labels
Ktree <- midpoint(Ktree)
write.nexus(Ktree,file=(paste(projectPath, "BAYESfiles/Ktree.nexus",sep="")),translate = TRUE)

#### Write files for each trait Bayestraits

BMAT <- read.csv(paste(projectPath, "Rfiles/boundMatrix.csv", sep="")) # matrix containing gene counts obtained with orthofinder
colnames(BMAT)<- gsub("\\.","-",colnames(BMAT))

SPEC <- read.csv(paste(projectPath, "SCOARYfiles/providedCSV.csv", sep=""))
colnames(SPEC)[1]<-"ID"

BMAT <- cbind("ID"=BMAT[,1],BMAT[,colnames(BMAT) %in% SPEC[,1]]) # verifies that the names in bound matrix are in the trait file

#Turn numeric columns into 0s and 1s
B1<-BMAT[,2:ncol(BMAT)]
B1<- as.data.frame(apply(B1, 2, function(x) ifelse(x > 1, 1, x))) # turn counts into presence/absence
BMAT<-cbind("ID"=BMAT[,1],B1)
BMAT <- BMAT[rowSums(BMAT[,2:(ncol(BMAT))])>1,] # eliminates rows where the trait is only in one genome
BMAT <- BMAT[rowSums(BMAT[,2:(ncol(BMAT))])<(ncol(BMAT)),]# eliminates rown where the trait is in ALL genomes


FORBT <- BMAT
FORBT <- as.data.frame(t(FORBT[,2:ncol(FORBT)])) # transpose the matrix
colnames(FORBT)<-BMAT$ID

FORBT$Name <- rownames(FORBT)
FORBT$Name <- gsub("\\.","-",FORBT$Name) #eliminate characters that interfere with downstream analysis
FORBT$Name <- gsub("\\-1","",FORBT$Name)
rownames(FORBT)<-FORBT$Name
FORBT <- subset(FORBT, select=-c(Name))
FORBT <- FORBT[!(row.names(FORBT) %in% ("ID")),]

TC <-  SPEC[,2] #METADATA COLUMN CONTAINING TRAIT TO COMPARE

write.csv(FORBT,paste(projectPath, "BAYESfiles/BayesGenerator.csv", sep=""),row.names = TRUE)
