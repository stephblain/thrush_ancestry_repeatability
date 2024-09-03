#module load GCC/11.2.0  OpenMPI/4.1.1 R/4.2.0

#R

library(tidyverse)

RUN_LOCATION="/scratch/user/sblain/get_AIMs/ref_panel/"
run_tag="ref"
sp1="inland"
sp2="coastal"

setwd(paste(RUN_LOCATION,run_tag,"_variantCalls",sep=""))

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
			mutate(ChromPos=paste(Chrom,Pos))%>%
			#Swainson's only - remove mtDNA
			filter(Chrom!="mtDNA")


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


#look at the number of 10kb windows included after filtering to different cutoffs

thinN<-10000 #number of sites to thin

for(i1 in rev(seq(0.3,0.9,by=0.1))){

panel1.i<-panel1.w%>%filter(sp2_sp1>i1)

sites1<-panel1.i%>%select(Chrom,Pos)
chrList<-unique(sites1$Chrom)


sites2<-data.frame()

for(x in chrList){
	chr1<-sites1%>%filter(Chrom==x)
	i=1
	while(i <= nrow(chr1)){
		j=chr1$Pos[i] #start position
		k=j+thinN #end position
		chr1<-chr1%>%filter(!(Pos>j&Pos<k)) #filter between
		i=i+1 #move to next row
	}
	sites2<-rbind(sites2,chr1)
}

print(paste("frequency level:",i1))
print(paste("number of 10kb windows:",nrow(sites2)))

}

##Outputs:
#
# "frequency level: 0.9"
# "number of 10kb windows: 15342"
# "frequency level: 0.8"
# "number of 10kb windows: 21643"
# "frequency level: 0.7"
# "number of 10kb windows: 34496"
# "frequency level: 0.6"
# "number of 10kb windows: 42987"
# "frequency level: 0.5"
# "number of 10kb windows: 53568"
# "frequency level: 0.4"
# "number of 10kb windows: 77805"
# "frequency level: 0.3"
# "number of 10kb windows: 87536"


#output 1: a list of AIMs that hit the frequency cutoff in ref panel 1
write.table(panel1.50%>%select(Chrom,Pos),
		file=paste("AIMs_f50",run_tag,"panel1",sep="_"),
		quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
write.table(panel1.70%>%select(Chrom,Pos),
		file=paste("AIMs_f70",run_tag,"panel1",sep="_"),
		quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
write.table(panel1.90%>%select(Chrom,Pos),
		file=paste("AIMs_f90",run_tag,"panel1",sep="_"),
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
	file=paste("parHMM_f50",run_tag,"panel2",sep="_"))
write.table(panel2.w%>%filter(ChromPos%in%panel1.70$ChromPos)%>%select(-ChromPos),
	file=paste("parHMM_f70",run_tag,"panel2",sep="_"))
write.table(panel2.w%>%filter(ChromPos%in%panel1.90$ChromPos)%>%select(-ChromPos),
	file=paste("parHMM_f90",run_tag,"panel2",sep="_"))



