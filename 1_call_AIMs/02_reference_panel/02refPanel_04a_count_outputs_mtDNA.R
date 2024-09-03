#module load GCC/11.2.0  OpenMPI/4.1.1 R/4.2.0

#R

library(tidyverse)

RUN_LOCATION="/scratch/user/sblain/get_AIMs/ref_panel/"
run_tag="ref"
add_tag="mtDNA"
sp1="inland"
sp2="coastal"

setwd(paste(RUN_LOCATION,run_tag,"_variantCalls",sep=""))

raw_aims<-read.table("/scratch/user/sblain/get_AIMs/define_AIMs/raw_AIMs_list",header=F)

#read in data
sp1p1<-read.table("ref_inland_panel1_AIMs_counts",header=F)
sp2p1<-read.table("ref_coastal_panel1_AIMs_counts",header=F)
sp1p2<-read.table("ref_inland_panel2_AIMs_counts",header=F)
sp2p2<-read.table("ref_coastal_panel2_AIMs_counts",header=F)

panel1<-rbind(sp1p1%>%mutate(species="species1"),
			sp2p1%>%mutate(species="species2"))%>%
			rename(Chrom=V1,Pos=V2,Total=V3,Alt=V4)%>%
			mutate(AltFreq=(Alt/Total))

#look at summary stats
panel1%>%group_by(species)%>%
	summarise(n_aims=n(),mean_alt=mean(Alt),
			mean_total=mean(Total),mean_AltFreq=mean(AltFreq,na.rm=T))
	
panel1.w<-panel1%>%
			select(Chrom,Pos,species,AltFreq)%>%
			pivot_wider(names_from=species,values_from=AltFreq)%>%
			mutate(sp2_sp1=species2-species1)%>%
			mutate(sp1_sp2=species1-species2)%>%
			mutate(ChromPos=paste(Chrom,Pos))
			
#check that the direction of sp2=Alt and sp1=Ref is correct
panel1.w%>%summarise(sp2_sp1_50=sum(sp2_sp1>0.5,na.rm=T),
						sp2_sp1_70=sum(sp2_sp1>0.7,na.rm=T),
						sp2_sp1_90=sum(sp2_sp1>0.9,na.rm=T),
						sp1_sp2_50=sum(sp1_sp2>0.5,na.rm=T),
						sp1_sp2_70=sum(sp1_sp2>0.7,na.rm=T),
						sp1_sp2_90=sum(sp1_sp2>0.9,na.rm=T))

panel1.50<-panel1.w%>%filter(sp2_sp1>0.5)


#output 1: a list of AIMs that hit the frequency cutoff in ref panel 1
write.table(panel1.50%>%select(Chrom,Pos),
		file=paste("AIMs_f50",run_tag,add_tag,"panel1",sep="_"),
		quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")


panel2<-rbind(sp1p2%>%mutate(species="species1"),
			sp2p2%>%mutate(species="species2"))%>%
			rename(Chrom=V1,Pos=V2,Total=V3,Alt=V4)%>%
			mutate(Ref=Total-Alt)
#look at summary stats
panel2%>%group_by(species)%>%
	summarise(n_aims=n(),mean_alt=mean(Alt),
			mean_total=mean(Total),mean_Ref=mean(Ref))

panel2.w<-panel2%>%select(-Total)%>%
		pivot_wider(names_from=species,values_from=c(Ref,Alt))%>%
		select(Chrom,Pos,Ref_species1,Alt_species1,Ref_species2,Alt_species2)%>%
		mutate(ChromPos=paste(Chrom,Pos))

#divide counts in half - haploid genome
panel2.w<-panel2.w%>%
		mutate(Ref_species1=ifelse(Chrom=="mtDNA",Ref_species1/2,Ref_species1),
			Ref_species2=ifelse(Chrom=="mtDNA",Ref_species2/2,Ref_species2),
			Alt_species1=ifelse(Chrom=="mtDNA",Alt_species1/2,Alt_species1),
			Alt_species2=ifelse(Chrom=="mtDNA",Alt_species2/2,Alt_species2))
		
#check for decimals		
panel2.w%>%filter(Chrom=="mtDNA"&
	(Ref_species1%%1==1|Ref_species1%%1==1|Alt_species2%%1==1|Alt_species2%%1==1))
		
colSums(panel2.w%>%filter(ChromPos%in%panel1.50$ChromPos)%>%
	select(Ref_species1,Alt_species1,Ref_species2,Alt_species2))/
	nrow(panel2.w%>%filter(ChromPos%in%panel1.50$ChromPos))

write.table(panel2.w%>%filter(ChromPos%in%panel1.50$ChromPos)%>%select(-ChromPos),
	file=paste("parHMM_f50",run_tag,add_tag,"panel2",sep="_"),
	quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")

#make new current aims file
current_aims<-raw_aims%>%
	mutate(V1V2=paste(V1,V2))%>%
	filter(V1V2%in%panel1.50$ChromPos)%>%
	select(-V1V2)

write.table(current_aims,file=paste("current_aims","f50",add_tag,sep="_"),
	quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
	



