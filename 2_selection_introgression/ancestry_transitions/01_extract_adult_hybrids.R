#module load GCC/12.2.0  OpenMPI/4.1.4 R/4.3.1

#Rscript --vanilla

#based on metadata, extract adult individuals from hybrid zones in Alaska, Pemberton, Hope, and Washington
#format ancestry states for input to ancestry transition estimation

library(tidyverse); library(data.table)

setwd("/scratch/user/sblain/genomic_clines/ancestryTracts")

aims<-fread("thrushHybrids.f50.240205.states",header=T)
meta<-fread("../AIMs_metadata_240215.csv",header=T)

aims<-aims%>%mutate(ChromPos=paste(CHROM,POS,sep="_"))%>%
  mutate(locusID=paste("locus",0:(nrow(aims)-1)))

hybMeta<-meta%>%
  filter(tag_type!="radio"&age_release!="HY"& #adults only
           !reference%in%c("bird_Z_plus","CF20H06","AF14H02")& #remove duplicate birds (bird_Z_plus=bird_Z,CF20H06=CF20H05,AF14H02=DF27H08)
           release_site%in%c("Alaska","Hope","Pemberton","Washington")&  #from hybrid zone
           !ref_taxon%in%c("inland","coastal"))%>% #not an inland or coastal ref
  filter(reference%in%colnames(aims))

runDay<-substr(gsub("-","",Sys.Date()),3,8)


###Make input for ancestry tracts

hybAims2<-aims%>%select(all_of(c("CHROM","POS","REF","ALT",hybMeta$reference)))
hybAims2out=paste("thrushHybrids.adultHybrids.f50",runDay,"states",sep=".")
write.table(hybAims2,file=hybAims2out)