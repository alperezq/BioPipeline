library(ape)
library(pheatmap)
library(stringr)
library(dplyr)
library(phangorn)
library(argparse)

args <- commandArgs(trailingOnly=TRUE)
projectPath <- args[1]

#set paths of files
KTreeFile <- paste(projectPath, "KSNP3files/kSNP3_results/tree.ML.tre", sep="")
BMATFile <- paste(projectPath, "Rfiles/boundMatrix.csv", sep="")
SPECFile <- paste(projectPath, "SCOARYfiles/providedCSV.csv", sep="")

KTreeFunc <- function()
{
  if(file.exists(KTreeFile))
  {
    KTree <- read.tree(KTreeFile)
    if(object.size(KTree)>0)
    {
      KTree$node.label<-NULL #delete node labels
      KTree <- midpoint(KTree)
      write.nexus(KTree,file=(paste(projectPath, "BAYESfiles/Ktree.nexus",sep="")),translate = TRUE)
    }else{write("tree.ML.tre file has insufficient data, unable to create Ktree.nexus.", file = err, append = TRUE)}
  }else{write("tree.ML.tre file was not found, unable to create Ktree.nexus.", file = err, append = TRUE)}
}

#### Write files for each trait Bayestraits
BoundMatFunc <- function()
{
  if(file.exists(BMATFile))
  {
    BMAT <- read.csv(BMATFile) # matrix containing gene counts obtained with orthofinder
    if(object.size(BMAT) > 0)
    {
      colnames(BMAT)<- gsub("\\.","-",colnames(BMAT))
      if(file.exists(SPECFile))
      {
        SPEC <- read.csv(SPECFile)
        if(object.size(SPEC) > 0)
        {
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

          TC <-  SPEC[,2] #METADATA COLUMN CONTAINING TRAIT TO COMPARE (not currently used)

          write.csv(FORBT,paste(projectPath, "BAYESfiles/BayesGenerator.csv", sep=""),row.names = TRUE)
        }else{write("Provided csv has insufficient data, unable to create BayesGenerator.csv.", file = err, append = TRUE)}
      }else{write("No provided csv found, unable to create BayesGenerator.csv.", file = err, append = TRUE)}
    }else{write("Bound Matrix has insufficient data, unable to create BayesGenerator.csv", file = err, append = TRUE)}
  }else{write("No Bound Matrix, unable to create BayesGenerator.csv.", file = err, append = TRUE)}
}

#Open error file for writing
errFile <- paste(projectPath,"Logging/RErrors.txt", sep = "")
err <- file(paste(projectPath,"Logging/RErrors.txt", sep = ""))

#call to above functions
KTreeFunc()
BoundMatFunc()

#Close error file
close(err)
