
Scoary<-read.csv("~/Desktop/BACTPIPES/Race_9b_03_07_2019_1244.results.csv")
Ids<-Scoary[Scoary$Bonferroni_p<0.05,1]


#############
A<-read.csv("/run/user/1000/gvfs/sftp:host=129.82.36.76/data/rex/BioPipeline/projects/TEST/Rfiles/RepeatNames.csv")
B<-read.delim("/run/user/1000/gvfs/sftp:host=129.82.36.76/data/rex/BioPipeline/projects/TEST/TRFfiles/trfParsed.txt",sep = " ",header = FALSE)

library(dplyr)

C<-left_join(A,B,by=c("Sequence"="V7"))

D<-C[C$Name %in% Ids,]

#write out D
################################## RVDs nucle

TALgroups<-read.csv("/run/user/1000/gvfs/sftp:host=129.82.36.76/data/rex/BioPipeline/projects/idunsTEST/DISTALfiles/Outputs/disTALOut.TALgroups.csv")
TalIds<-TALgroups$TAL[paste("TALGroup_",TALgroups$Group,sep="") %in% Ids]

A<-read.delim("/run/user/1000/gvfs/sftp:host=129.82.36.76/data/Alvaro/Concatenated_createdbyalvaro.csv",sep="\t")
B <- A[A$ID %in% TalIds,]

library(ape)

C<-read.FASTA("/run/user/1000/gvfs/sftp:host=129.82.36.76/data/rex/BioPipeline/projects/TEST/DISTALfiles/rvdCombo.FASTA")
C<-C[names(C) %in% TalIds]

write.FASTA(C,"~/Desktop/Blabla.fasta")

################################## Orthofinder
library(reshape2)
library(stringr)
A<-read.delim("/run/user/1000/gvfs/sftp:host=129.82.36.76/data/rex/BioPipeline/projects/TEST/ORTHOfiles/Orthogroups.txt",sep=" ",header = FALSE)

B <- melt(A,id.vars = "V1")
B <- B[,c(1,3)]
B <- B[B$value != "",]
B$V1 <- gsub(":", "",B$V1)



## note to create the cocatenated file, add the genome naime to each header in the proteins file

#>/data/Alvaro/AllFAAs_withnames.fa
#for i in *.faa
#do
#NAM=$(echo $i| sed 's/.faa//g')                   
#sed "s/>/>$NAM /g" $i >>/data/Alvaro/AllFAAs_withnames.fa
#done


D<-read.FASTA("/run/user/1000/gvfs/sftp:host=129.82.36.76/data/Alvaro/AllFAAs_withnames.fa")
fastaIds<-data.frame("Short_ID"=str_split_fixed(names(D)," ",3)[,2],"fastaID"=names(D))
E<-left_join(fastaIds,B,by=c("Short_ID"="value"))
E$New_ID<-paste(E$V1,E$fastaID,sep=":")

names(D)<-E$New_ID

H <- sort(E$New_ID[E$V1 %in% Ids])

I <- D[names(D) %in% H]

J <- I[order(names(I))]
write.FASTA(J,"~/Desktop/SIGOrthogroups.fa")
  