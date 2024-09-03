#!/bin/bash
      
##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=make_consensus    	#Set the job name
#SBATCH --time=4:00:00          	#Set the wall clock limit
#SBATCH --ntasks=1               	#Request 1 task
#SBATCH --mem=8G              		#Request xGB per task
#SBATCH --output=error_05_consensus.%j	#Send stdout/err to filename
#SBATCH --gres=gpu:1             	#Request 1 GPU per node can be 1 or 2
#SBATCH --partition=gpu          	#Request the GPU partition/queue


#This filters the called variants then uses that to make a consensus sequence with the reference
# so that there is a pseudoref that is collinear with the inland ref


cd /scratch/user/sblain/get_AIMs/define_AIMs

module load GCC/11.2.0
module load BCFtools/1.14

# perform the filtering with bcftools

bcftools view --exclude-types indels coastal.bamOut.sorted.unique.vcf.gz | \ #no indels
 bcftools view --genotype hom | \ #homozygous only (will also filter to only biallelic for this dataset)
 bcftools view -e 'QUAL < 30 || DP < 7 || DP > 60 || MQ < 40 ' \ #other assorted thresholds
 -o coastal.bamOut.sorted.unique.filtered.vcf.gz




