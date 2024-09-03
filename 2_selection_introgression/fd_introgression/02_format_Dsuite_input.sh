#!/bin/bash

#SBATCH --job-name=format_dsuite_input
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --output=error_format.%j

module purge
module load GCC/11.2.0 BCFtools/1.14 VCFtools/0.1.16

vcf=/scratch/user/sblain/thrush_hybrids/d_stats/variantCalling/thrush.d_stats.filtered.vcf.gz
prefix=thrushDstats

#remove Bicknell's thrush - keep Hermit as the outgroup
bcftools view -s ^KSW3633 -m2 -M2 -v snps $vcf | \
        bcftools filter -e 'F_MISSING > 0.25 || MAF <= 0.05' -O z -o $prefix.filtered.vcf.gz

tabix $prefix.filtered.vcf.gz