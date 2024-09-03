#!/usr/bin/Rscript --vanilla

#modules to run in on cluster: module load GCC/12.2.0 OpenMPI/4.1.4 R/4.3.1
#Run in cluster as: Rscript --vanilla run_bgchm_thrushes_withZ.R population_name
#population_name should be Alaska, Pemberton, Hope, or Washington

#run once to install (takes a while):
#library(devtools)
##not enough space in home directory for me to install there - need to install in scratch
#library(withr) #for with_libpaths(), needed to specify directory
#with_libpaths(new = "/scratch/user/sblain/tools", install_github('zgompert/bgc-hm'))

# #To test installation, run:
# ## load the data set
# data(genotypes)
# ## this includes three objects, GenHybrids, GenP0, and GenP1
# ## estimate parental allele frequencies, uses default HMC settings
# p_out<-est_p(G0=GenP0,G1=GenP1,model="genotype",ploidy="diploid")
# 
# ## estimate hybrid indexes, uses default HMC settings
# ## and uses point estimates (posterior medians) of allele frequencies
# h_out<-est_hi(Gx=GenHybrids,p0=p_out$p0[,1],p1=p_out$p1[,1],model="genotype",ploidy="diploid")

#Note: Washington had trouble initializing, had to restart script 3 times

library("bgchm", lib.loc = '/scratch/user/sblain/tools/') #loads from folder in this directory called bgchm
library(tidyverse)

#run on command line, in $SCRATCH/tools: git clone https://github.com/zgompert/bgc-hm.git
#run in R: devtools::install("/scratch/user/sblain/tools/bgc-hm")
setwd("/scratch/user/sblain/genomic_clines/bgchm_clines")

#read in aims metadata to match individual IDs to populations
meta<-read_csv("AIMs_metadata_240227.csv")


aims<-read_table("thrushHybrids.f50.240205.LD.probs")
aimsStates<-read_table("thrushHybrids.f50.240205.LD.states")
Zaims<-read_table("thrushHybridsZchrom.f50.240307.LD.probs")
ZaimsStates<-read_table("thrushHybridsZchrom.f50.240307.LD.states")


args=commandArgs(trailingOnly=TRUE)
pop1=args[1] #pop1 = "Washington"
print(paste("focal population:",pop1))

#filter metadata to only include hybrids to fit clines to
pop1Meta<-meta%>%
  filter(tag_type!="radio"&age_release!="HY"& #adults only
           !reference%in%c("bird_Z_plus","CF20H06","AF14H02")& #remove duplicate birds (bird_Z_plus=bird_Z,CF20H06=CF20H05,AF14H02=DF27H08)
           release_site==pop1&  #from focal population
           !ref_taxon%in%c("inland","coastal"))%>% #not an inland or coastal ref
  filter(reference%in%colnames(Zaims))

#get list of IDs for parent 0
par0IDs<-meta%>%filter(ref_taxon=="inland"&reference%in%colnames(Zaims))%>%pull(reference)
allAimsPar0<-rbind(aimsStates%>%mutate(ChromPos=paste(CHROM,POS,sep="_"))%>%
                     select(all_of(c("ChromPos",par0IDs))),
                   ZaimsStates%>%mutate(ChromPos=paste(CHROM,POS,sep="_"))%>%
                     select(all_of(c("ChromPos",par0IDs))))

#allele frequency for inlands at each marker
#don't allow zero (i.e. fixed SNP) - replace with 0.001
#fixes computational problem, and likely more biologically realistic anyway
par0Aims<-allAimsPar0%>%pivot_longer(cols=-ChromPos)%>%
  group_by(ChromPos)%>%summarise(genotype=mean(value/2,na.rm=T))%>%
  mutate(genotype=if_else(genotype<0.001,0.001,genotype))

#get list of IDs for parent 1
par1IDs<-meta%>%filter(ref_taxon=="coastal"&reference%in%colnames(Zaims))%>%pull(reference)
allAimsPar1<-rbind(aimsStates%>%mutate(ChromPos=paste(CHROM,POS,sep="_"))%>%
                     select(all_of(c("ChromPos",par1IDs))),
                   ZaimsStates%>%mutate(ChromPos=paste(CHROM,POS,sep="_"))%>%
                     select(all_of(c("ChromPos",par1IDs))))

#allele frequency for coastals
#don't allow fixed SNP - replace with 0.999
par1Aims<-allAimsPar1%>%pivot_longer(cols=-ChromPos)%>%
  group_by(ChromPos)%>%summarise(genotype=mean(value/2,na.rm=T))%>%
  mutate(genotype=if_else(genotype>0.999,0.999,genotype))


aims<-aims%>%mutate(ChromPos=paste(CHROM,POS,sep="_"))
Zaims<-Zaims%>%mutate(ChromPos=paste(CHROM,POS,sep="_"))

allAims=rbind(aims%>%select(all_of(c("ChromPos",pop1Meta$reference))),
              Zaims%>%select(all_of(c("ChromPos",pop1Meta$reference))))

#extract geno probabilities and convert to numeric
pop1Aims<-allAims%>%
  pivot_longer(cols=-ChromPos,names_to="reference")%>%
  mutate(g1pr=as.numeric(str_split_i(value,",",1)),
         g2pr=as.numeric(str_split_i(value,",",2)),
         g3pr=as.numeric(str_split_i(value,",",3)))

#get a matrix of posterior probabilities associated with each ancestry state
g1mat<-pop1Aims%>%select(ChromPos,reference,g1pr)%>%
  pivot_wider(names_from = ChromPos,values_from = g1pr)%>%
  column_to_rownames(var="reference")%>%
  as.matrix()
g2mat<-pop1Aims%>%select(ChromPos,reference,g2pr)%>%
  pivot_wider(names_from = ChromPos,values_from = g2pr)%>%
  column_to_rownames(var="reference")%>%
  as.matrix()
g3mat<-pop1Aims%>%select(ChromPos,reference,g3pr)%>%
  pivot_wider(names_from = ChromPos,values_from = g3pr)%>%
  column_to_rownames(var="reference")%>%
  as.matrix()
genoList<-list(g1mat,g2mat,g3mat)

print(paste("number of hybrid individuals:",dim(g1mat)[1]))
print(paste("number of loci:",dim(g1mat)[2]))

#reorder parental loci to match hybrids
par0Aims<-par0Aims[match(colnames(g1mat),par0Aims$ChromPos),]
par1Aims<-par1Aims[match(colnames(g1mat),par1Aims$ChromPos),]

print(paste("mean parent 0 genotype:", round(mean(par0Aims$genotype),4)))
print(paste("mean parent 1 genotype:", round(mean(par1Aims$genotype),4)))

#check that loci (colnames) and individuals (rownames) are in the same order in input files
#if not - don't fit model, print error
if(all(colnames(g1mat)==par0Aims$ChromPos)&
   all(colnames(g1mat)==par1Aims$ChromPos)&
   all(colnames(g1mat)==colnames(g2mat))&
   all(colnames(g1mat)==colnames(g3mat))&
   all(rownames(g1mat)==pop1Meta$reference)&
   all(rownames(g1mat)==rownames(g2mat))&
   all(rownames(g3mat)==rownames(g3mat))){
  
  #fit hierarchical genotype likelihood model
  fit1<-est_genocl(Gx=genoList,p0=par0Aims$genotype,p1=par1Aims$genotype,
                   H=pop1Meta$aims_ancestry,
                   model="glik",ploidy="diploid",hier=TRUE,n_iters=10000)
  
  #format output cline centres
  centreOut<-data.frame(fit1$center)%>%
    rename(median="X50.",CIlwr="X5.",CIupr="X95.")%>%
    mutate(locus=colnames(g1mat),population=pop1,param="centre")
  
  #format output cline gradients
  gradientOut<-data.frame(fit1$gradient)%>%
    rename(median="X50.",CIlwr="X5.",CIupr="X95.")%>%
    mutate(locus=colnames(g1mat),population=pop1,param="gradient")
  
  pop1Out<-rbind(centreOut,gradientOut)
  
  print(fit1$gencline_hmc)
  
  
  write.table(pop1Out,file=paste("bgchm_clines_out_withZ",pop1,"tsv",sep="."),
              quote=F,sep="\t",col.names=T,row.names=F)
  
}else{print(paste("ERROR IN INPUT FILES: Check dataframe matching for",pop1))}

