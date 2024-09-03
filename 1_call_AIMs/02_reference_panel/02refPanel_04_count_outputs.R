#module load GCC/11.2.0  OpenMPI/4.1.1 R/4.2.0

#R

library(tidyverse)

RUN_LOCATION="/scratch/user/sblain/get_AIMs/ref_panel/"
run_tag="ref"
sp1="inland"
sp2="coastal"

setwd(paste(RUN_LOCATION,run_tag,"_variantCalls",sep=""))

#read in data
sp1<-read.table("ref_inland_AIMs_counts",header=F)
sp2<-read.table("ref_coastal_AIMs_counts",header=F)

panel1<-rbind(sp1%>%mutate(species="species1"),
			sp2%>%mutate(species="species2"))%>%
			rename(Chrom=V1,Pos=V2,Total=V3,Alt=V4,RefAllele=V5,AltAllele=V6)%>%
			mutate(AltFreq=(Alt/Total))

#look at summary stats
panel1%>%group_by(species)%>%
	summarise(n_aims=n(),mean_alt=mean(Alt),
			mean_total=mean(Total),mean_AltFreq=mean(AltFreq,na.rm=T))
	
panel1.w<-panel1%>%
			select(Chrom,Pos,species,AltFreq,RefAllele,AltAllele)%>%
			pivot_wider(names_from=species,values_from=AltFreq)%>%
			mutate(sp2_sp1=species2-species1)%>%
			mutate(sp1_sp2=species1-species2)%>%
			mutate(ChromPos=paste(Chrom,Pos))#%>%
			#Swainson's only - remove mtDNA
			#filter(Chrom!="mtDNA")
			
#check that the direction of sp2=Alt and sp1=Ref is correct
panel1.w%>%summarise(sp2_sp1_50=sum(sp2_sp1>0.5,na.rm=T),
						sp2_sp1_70=sum(sp2_sp1>0.7,na.rm=T),
						sp2_sp1_90=sum(sp2_sp1>0.9,na.rm=T),
						sp1_sp2_50=sum(sp1_sp2>0.5,na.rm=T),
						sp1_sp2_70=sum(sp1_sp2>0.7,na.rm=T),
						sp1_sp2_90=sum(sp1_sp2>0.9,na.rm=T))

panel1.50<-panel1.w%>%filter(sp2_sp1>0.5)
panel1.70<-panel1.w%>%filter(sp2_sp1>0.7)
panel1.90<-panel1.w%>%filter(sp2_sp1>0.9)


#output 1: a list of AIMs that hit the frequency cutoff in ref panel 1, with underscores in scaffold names
write.table(panel1.50%>%select(Chrom,Pos),
		file="AIMs_f50.pos",
		quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
write.table(panel1.70%>%select(Chrom,Pos),
		file="AIMs_f70.pos",
		quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
write.table(panel1.90%>%select(Chrom,Pos),
		file="AIMs_f90.pos",
		quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")


#4 tab-delimited columns: chromosome, position, ref allele, alt allele
#in same folder, have aims.mod file with tab-delimited: chromosome_position, position, ref allele, alt allele
# and aims.mod.bed file with tab-delimited: chromosome, position, position, ref allele, alt allele
panel1.50.aims<-panel1.50%>%mutate(Chrom=gsub("_","-",Chrom))
write.table(panel1.50.aims%>%select(Chrom,Pos,RefAllele,AltAllele),
		file="AIMs_f50",
		quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")

panel1.50.mod<-panel1.50.aims%>%mutate(ChromPos=paste(Chrom,Pos,sep="_"))%>%
	select(ChromPos,Pos,RefAllele,AltAllele)
write.table(panel1.50.mod,
		file="AIMs_f50.mod",
		quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")

panel1.50.bed<-panel1.50.aims%>%mutate(Pos1=Pos)%>%select(Chrom,Pos,Pos1,RefAllele,AltAllele)
write.table(panel1.50.bed,
		file="AIMs_f50.mod.bed",
		quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")





panel2<-rbind(sp1%>%mutate(species="species1"),
			sp2%>%mutate(species="species2"))%>%
			rename(Chrom=V1,Pos=V2,Total=V3,Alt=V4,RefAllele=V5,AltAllele=V6)%>%
			mutate(Ref=Total-Alt)
#look at summary stats
panel2%>%group_by(species)%>%
	summarise(n_aims=n(),mean_alt=mean(Alt),
			mean_total=mean(Total),mean_Ref=mean(Ref))

panel2.w<-panel2%>%select(-Total)%>%
		pivot_wider(names_from=species,values_from=c(Ref,Alt))%>%
		select(Chrom,Pos,Ref_species1,Alt_species1,Ref_species2,Alt_species2)%>%
		mutate(ChromPos=paste(Chrom,Pos))%>% #chrompos has underscores in scaffold names
		mutate(Chrom=gsub("_","-",Chrom)) #no underscores in scaffold names
		
colSums(panel2.w%>%filter(ChromPos%in%panel1.50$ChromPos)%>%
	select(Ref_species1,Alt_species1,Ref_species2,Alt_species2))/
	nrow(panel2.w%>%filter(ChromPos%in%panel1.50$ChromPos))
colSums(panel2.w%>%filter(ChromPos%in%panel1.70$ChromPos)%>%
	select(Ref_species1,Alt_species1,Ref_species2,Alt_species2))/
	nrow(panel2.w%>%filter(ChromPos%in%panel1.70$ChromPos))
colSums(panel2.w%>%filter(ChromPos%in%panel1.90$ChromPos)%>%
	select(Ref_species1,Alt_species1,Ref_species2,Alt_species2))/
	nrow(panel2.w%>%filter(ChromPos%in%panel1.90$ChromPos))

write.table(panel2.w%>%filter(ChromPos%in%panel1.50$ChromPos)%>%select(-ChromPos),
	file="parHMM_f50",
	quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
write.table(panel2.w%>%filter(ChromPos%in%panel1.70$ChromPos)%>%select(-ChromPos),
	file="parHMM_f70",
	quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
write.table(panel2.w%>%filter(ChromPos%in%panel1.90$ChromPos)%>%select(-ChromPos),
	file="parHMM_f90",
	quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")



