#!/bin/bash

#SBATCH --job-name=ref_counts
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --output=error_ref_counts.%j

module load GCC/11.2.0 BCFtools/1.14


RUN_LOCATION="/scratch/user/sblain/get_AIMs/ref_panel/"
run_tag="ref"
IND_LIST="/scratch/user/sblain/get_AIMs/ref_panel/ind_list_ref"

IND_LIST_sp1="/scratch/user/sblain/get_AIMs/ref_panel/inds_ref_inland"
IND_LIST_sp2="/scratch/user/sblain/get_AIMs/ref_panel/inds_ref_coastal"

sp1="inland"
sp2="coastal"



cd $RUN_LOCATION/"$run_tag"_variantCalls

#concatenate all AIMs files and only keep biallelic snps
bcftools concat *AIMs.vcf.gz | \
	bcftools view -m2 -M2 -v snps \
		-o all_scaffolds."$run_tag".AIMs.vcf.gz

tabix all_scaffolds."$run_tag".AIMs.vcf.gz

echo "AIMs vcf concatenated"

#bcftools select coastal individuals and get their allele counts
bcftools view -S $IND_LIST_sp1 all_scaffolds."$run_tag".AIMs.vcf.gz | \
bcftools query -f '%CHROM %POS %AN %AC %REF %ALT\n' > "$run_tag"_"$sp1"_AIMs_counts

bcftools view -S $IND_LIST_sp2 all_scaffolds."$run_tag".AIMs.vcf.gz | \
bcftools query -f '%CHROM %POS %AN %AC %REF %ALT\n' > "$run_tag"_"$sp2"_AIMs_counts


echo "ref allele counts estimated"


