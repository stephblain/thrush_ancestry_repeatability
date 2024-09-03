#use inland genome annotation to estimate the number of genes per 100kb window

library(tidyverse)
library(genomation) # I was not able to load all these libraries, so i found other programs to load gff files online
library(GenomicRanges)  # I was not able to load all these libraries, so i found other programs to load gff files online

setwd("/scratch/user/sblain/thrush_hybrids/d_stats/")

gff <- gffToGRanges("/scratch/user/sblain/genomic_clines/get_genes/GCF_009819885.2_bCatUst1.pri.v2_genomic_wfxnl.gff")
gff_genes <- subset(gff,type=="gene")

chroms<-read.table("/scratch/user/sblain/get_AIMs/define_AIMs/inland_ref/inland_genome.fasta.fai")
chroms<-chroms%>%select(V1,V2)%>%
  rename("V1"="scaffold","V2"="scaffoldLength")

scafOrder<-read.csv("/scratch/user/sblain/genomic_clines/scaffold_order.csv")
scafOrder<-scafOrder%>%as_tibble()%>%
  rename("RefSeq.Accn"="seqnames","Sequence.Name"="scaffold")%>%
  select(scaffold,seqnames)

chroms<-chroms%>%left_join(scafOrder)%>%filter(!is.na(seqnames))
gff_genes
gene_densities<-data.frame()

for(seq.i in unique(chroms$seqnames)){
  scafLength.i=chroms%>%filter(seqnames==seq.i)%>%pull(scaffoldLength)
  start.i=1
  end.i=100000
  while(start.i<scafLength.i){

    q=GRanges(seqnames=seq.i,
              ranges=IRanges(start=start.i, end=end.i))
    gene_densities<-rbind(gene_densities,data.frame(seqnames=seq.i,startPos=start.i,endPos=end.i,
                                                    n_genes=length(findOverlaps(gff_genes,q))))

    start.i=start.i+100000
    end.i=end.i+100000
  }
}

options(scipen = 20)

#were all chroms processed?
all(unique(gene_densities$seqnames)==unique(chroms$seqnames))

gene_densities<-left_join(gene_densities,chroms)
gene_densities<-gene_densities%>%mutate(endPos=if_else(scaffoldLength<endPos,scaffoldLength,endPos))%>%
  mutate(gene_density=n_genes/(endPos-startPos))

write.csv(gene_densities,file="gene_densities_100kb.csv",row.names=F)