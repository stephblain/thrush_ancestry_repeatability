#!/bin/bash

module purge
module load GCC/11.2.0 BCFtools/1.14 VCFtools/0.1.16

vcf=/scratch/user/sblain/thrush_hybrids/d_stats/variantCalling/thrush.d_stats.filtered.vcf.gz
prefix=thrushDstats

$SCRATCH/tools/Dsuite/Build/Dsuite Dinvestigate -w 100,100 $prefix.filtered.vcf.gz $prefix.species.tsv $prefix.trios.tsv