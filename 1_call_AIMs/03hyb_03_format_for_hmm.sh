
#################
#Run in bash
#################

###format hmm input file

IND_LIST="/scratch/user/sblain/get_AIMs/hmm_AIMs/inds_lowCov_2023" #list of individuals in order
Curr_AIMS="current_aims_f50.cMs" #with genetic distance in cMs in third column
RUN_LOCATION="/scratch/user/sblain/get_AIMs/hmm_AIMs/lowCov/lowCov_hmm_f50" #location of input and output files
RAW_COUNTS="/scratch/user/sblain/get_AIMs/hmm_AIMs/lowCov/lowCov_variantCalls/lowCov.hmmFormat_all"
PAR_FILE="parHMM_f50_ref_panel2"
run_tag="lowCov" #use the same run tag across all scripts
f_tag="f50" #frequency cutoff

cd $RUN_LOCATION

sed "s/\t/_/" $Curr_AIMS > $Curr_AIMS.mod
cut -f 1 $Curr_AIMS.mod > "$Curr_AIMS".mod_1


#note - lowCov.hmmFormat_all includes mtDNA

#grep lines that match pattern 
#w=word matching; f=get pattern from file; F=compare strings (not regular expressions)
grep -wFf "$Curr_AIMS".mod_1 $RAW_COUNTS > "$run_tag"."$f_tag".hmmFormat_all

sed "s/\t/_/" $PAR_FILE > $PAR_FILE.mod
grep -wFf "$Curr_AIMS".mod_1 $PAR_FILE.mod > $PAR_FILE.$run_tag.$f_tag.mod
sed "s/_/\t/" $PAR_FILE.$run_tag.$f_tag.mod > $PAR_FILE.$run_tag.$f_tag



#sort by chrom and position before joinging
#sorts by col 1 as character then 2 as numeric
sort -k 1,1 -k 2,2n $PAR_FILE.$run_tag.$f_tag > $PAR_FILE.$run_tag.$f_tag.sorted
sed "s/_/\t/" "$run_tag"."$f_tag".hmmFormat_all | sort -k 1,1 -k 2,2n > "$run_tag"."$f_tag".hmmFormat_all.sorted
sort -k 1,1 -k 2,2n $Curr_AIMS > $Curr_AIMS.sorted
paste $PAR_FILE.$run_tag.$f_tag.sorted $Curr_AIMS.sorted "$run_tag"."$f_tag".hmmFormat_all.sorted | sed "s/ /\t/g" | cut -f 1-6,11,14- > "$f_tag".hmm.counts




###Get sample list and format for hmm


ls *hmmFormat > ../"run_tag"_"f_tag"_fileList

module load GCC/11.2.0  OpenMPI/4.1.1 R/4.2.0

cd $RUN_LOCATION

R

#################
#Run in R
#################

library(tidyverse)
run_tag<-"lowCov"
tab1<-read.table(paste(run_tag,"fileList",sep="_"),header=F)
tab1<-tab1%>%mutate(V1=gsub(".hmmFormat","",V1),
	V2=2)
write.table(tab1,file=paste(run_tag,"samples",sep="_"),quote=FALSE, 
            col.names=FALSE,row.names=FALSE,sep="\t")
#

##########pre-hmm filtering

tab2<-read.table("f50.hmm.counts",header=F)
#tab2[tab2==0]<-NA
#tab2Rows<-rowSums(is.na(tab2[,8:ncol(tab2)]))
tab2Rows<-rowSums(tab2[,8:ncol(tab2)])

tab2Out<-tab2[tab2Rows>100,] #set a very reasonable cutoff
write.table(tab2Out,file="f50.hmm.counts.filtered",quote=FALSE, 
            col.names=FALSE,row.names=FALSE,sep="\t")