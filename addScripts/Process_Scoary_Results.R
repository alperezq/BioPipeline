library(ape)
library(pheatmap)
library(phangorn)
library(dendextend)

Scoary <- read.csv("~/Desktop/Race_9b_03_07_2019_1244.results.csv")
Sig_Ids<-Scoary$V1[Scoary$Bonferroni_p < 0.05]

boundmat <- read.csv("~/Desktop/boundMatrix_3.csv")
Sigmat<-boundmat[boundmat$V1 %in% Sig_Ids,2:ncol(boundmat)]
rownames(Sigmat)<-boundmat[boundmat$V1 %in% Sig_Ids,1]
pheatmap(Sigmat)

Ktree <- read.tree("/run/user/1000/gvfs/sftp:host=129.82.36.76/data/rex/pipeScripts/TEST/KSNP3files/kSNP3_results/tree.ML.tre")
Ktree$tip.label<-gsub("-","\\.",Ktree$tip.label)

Ktree<-midpoint(Ktree)
Ktree<-chronos(Ktree)
#Kclus<-as.hclust(Ktree)
mydend <- as.dendrogram(Ktree)
clade_order <- order.dendrogram(as.dendrogram(mydend))
clade_name <- labels(mydend)
clade_position <- data.frame(clade_name,
                             clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name,
                   colnames(Sigmat))
mat2 <- Sigmat[,new_order]
Kclus<-as.hclust(mydend)

Trait<-read.csv("~/Desktop/Trait.csv")
tannot<-data.frame("Trait"=Trait[,2])
rownames(tannot)<-Trait[,1]

pheatmap(mat2,cluster_cols = Kclus,annotation_col = tannot,treeheight_col = 200)
