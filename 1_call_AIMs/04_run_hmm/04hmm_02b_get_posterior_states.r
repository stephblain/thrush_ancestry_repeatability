#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
run_tagR = as.character(args[1])

#print state predicted by hmm with a posterior probability of 0.9 or greater

IND_LIST<-paste("hmm_inds",run_tagR,sep="_")

library(tidyverse)

ind_list<-read.table(IND_LIST,header=F)


for(IND in ind_list$V1){
  probs<-read.table(paste(IND,"posterior",sep="."),header=T)
  probs.l<-probs%>%pivot_longer(cols=c("X2.0","X1.1","X0.2"),
                                names_to="state.id",values_to="probs")%>%
    mutate(state=case_when(probs>0.9~state.id))%>%drop_na()%>%
    mutate(state=case_when(state=="X2.0"~0,state=="X1.1"~1,state=="X0.2"~2))
  probs.out<-left_join(probs,probs.l)%>%select(state)
  colnames(probs.out)<-IND
  write.table(probs.out,paste(IND,"posterior_state",sep="."),quote=FALSE,
              col.names=TRUE,row.names=FALSE,sep="\t")}
#