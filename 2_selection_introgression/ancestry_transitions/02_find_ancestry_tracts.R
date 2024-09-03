#modules to load: module load GCC/12.2.0  OpenMPI/4.1.4 R/4.3.1
#run as: Rscript --vanilla ancestry_tracts.R

library(tidyverse)
theme_set(theme_classic())

#read in ancestry states
#first few columns - "CHROM" "POS" "REF" "ALT" - describe each locus, subsequent columns are per individual
#ancestry states are 0, 1, 2, or NA
aims<-read.table("thrushHybrids.adultHybrids.f50.240222.states",header=T)
inds<-colnames(aims%>%select(-CHROM,-POS,-REF,-ALT))

chroms<-names(table(aims$CHROM))[table(aims$CHROM)>20]

#start the start and stop positions of ancestry tracts
tracts<-data.frame()

for(scaf in chroms){
  
  print(scaf)
  
  for(ind in inds){
    
    ind1<-aims%>%filter(CHROM==scaf)%>%select(POS,all_of(ind))%>%
      drop_na()
    
    if(nrow(ind1)>2){ #require at least three aims on the chromosome
    
    count1<-0
    tracts1<-data.frame(CHROM=scaf,reference=ind,start=ind1$POS[1],stop=NA,ancestry=ind1[1,2])
    
    
    for(i in 2:nrow(ind1)){
      if(ind1[i,2]!=ind1[i-1,2]){
        count1=count1+1
        tracts1[count1,"stop"]<-ind1[i,"POS"]
        tracts1<-rbind(tracts1,data.frame(CHROM=scaf,reference=ind,
                                          start=ind1[i,"POS"],stop=NA,ancestry=ind1[i,2]))
      } } 
    
    tracts1[nrow(tracts1),"stop"]<-ind1[nrow(ind1),"POS"]
    tracts<-rbind(tracts,tracts1)
  } } }

chrom_sum<-aims%>%group_by(CHROM)%>%filter(CHROM%in%chroms)%>%
  summarise(max_pos=max(POS),n_aims=n())


tracts<-tracts%>%mutate(tractLength=stop-start)%>%
  filter(tractLength>10000)%>% #remove tracts shorter than 10kb
  left_join(chrom_sum)%>%
  mutate(length_per_aims=tractLength/n_aims)%>%
  mutate(heterozygous=if_else(ancestry==1,1,0))
tractStarts<-tracts%>%group_by(CHROM,start)%>%
  summarise(start_count=n())%>%
  mutate(start_Mb=start/1000000,yVar=0)

write.csv(tracts,"tracts_ThrushHybridsAdults.csv",row.names=F)

#make some summary plots

gg1<-ggplot(tractStarts,
           aes(x=start_Mb,y=start_count))+
  geom_point(fill="grey10")+
  geom_point(aes(x=start_Mb,y=yVar),colour="firebrick",shape="|")+
  facet_grid(cols=vars(CHROM),space="free",scales="free")

ggsave("startPositionCounts_ThrushHybridsAdults.png",gg1,width=24,height=4)

gg2<-ggplot(tracts,aes(x=start/1000000,fill=CHROM))+
  geom_histogram()+
  scale_fill_manual(values=rep(c("grey10","grey60"),4))+
  theme(legend.position="none")+
  xlab("Position (Mb)")+
  facet_grid(cols=vars(CHROM),space="free",scales="free")

ggsave("startPositionHistogram_ThrushHybridsAdults.png",gg2,width=24,height=4)

aims2<-aims%>%filter(CHROM%in%chroms)

#estimate the length of each ancestry tract
tractLengths<-data.frame()

for(i in 1:nrow(aims2)){
  
  chr1<-aims2[i,"CHROM"]
  pos1<-aims2[i,"POS"]
  
  tractLengths<-rbind(tractLengths,
                      tracts%>%filter(CHROM==chr1)%>%
                        filter(start<pos1&stop>pos1)%>%
                        summarise(tractMean=mean(tractLength),
                                  ancestryMean=mean(ancestry),
                                  heterozygosity=sum(heterozygous)/n())%>%
                        mutate(CHROM=chr1,POS=pos1))
}


write.csv(tractLengths,"tractLengths_ThrushHybridsAdults.csv",row.names=F)

#make some summary plots

gg3<-ggplot(tractLengths,aes(x=POS/1000000,y=tractMean/1000000))+
  geom_point(aes(colour=CHROM))+xlab("Position (Mb)")+ylab("Mean tract length")+
  scale_colour_manual(values=rep(c("grey10","grey60"),50))+
  theme(legend.position="none")+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  facet_grid(cols=vars(CHROM),space="free",scales="free")

ggsave("meanTractLengths_ThrushHybridsAdults.png",gg3,width=24,height=4)


gg4<-ggplot(tractLengths,aes(x=POS/1000000,y=heterozygosity))+
  ylim(0,1)+
  geom_point(aes(colour=CHROM))+xlab("Position (Mb)")+ylab("Heterozygosity")+
  scale_colour_manual(values=rep(c("grey10","grey60"),50))+
  theme(legend.position="none")+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  facet_grid(cols=vars(CHROM),space="free",scales="free")

ggsave("meanHeterozygosity_ThrushHybridsAdults.png",gg4,width=24,height=4)

