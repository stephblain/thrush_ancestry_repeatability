#!/bin/bash
      
##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=06_seqtk     	#Set the job name
#SBATCH --time=6:00:00          	#Set the wall clock limit to x hrs
#SBATCH --ntasks=1               	#Request 1 task
#SBATCH --mem=128G             		#Request 12GB per task
#SBATCH --output=error_06_seqtk.%j	#Send stdout/err to "Example4Out.[jobID]"
#SBATCH --gres=gpu:1             	#Request 1 GPU per node can be 1 or 2
#SBATCH --partition=gpu          	#Request the GPU partition/queue

cd /scratch/user/sblain/get_AIMs/define_AIMs

module load GCC/10.2.0
module load seqtk/1.3

seqtk mutfa inland_genome.fasta \
  coastal.bamOut.sorted.unique.vcf.cov-corrected.insnp > coastal.pseudoref.fasta


#NOTE: If you run this, follow it with picard Normalize to make characters per line the same as in the reference
#ALSO: If you run this, request a decent amount of memory. It doesn't take long to run, with this, but will crash with less

module load picard/2.25.1-Java-11

#normalize seqtk output to have the same line lengths as inland reference

java -jar $EBROOTPICARD/picard.jar NormalizeFasta I=coastal.pseudoref.fasta O=coastal.pseudoref.normalized.fasta LINE_LENGTH=50



#Alternative: make a consensus sequence
#this does the same thing as seqtk mutfa but uses less memory and requires less formatting
cat inland_genome.fasta | bcftools consensus coastal.bamOut.sorted.unique.filtered.vcf.gz > coastal_consensus.fasta


#identify_AIMs_two_genomes_v2.pl is a script downloaded from the Schumer lab github
#finds and prints sites that are different between the collinear genomes
perl scripts/Schumer_scripts/identify_AIMs_two_genomes_v2.pl inland_genome.fasta coastal.pseudoref.normalized.fasta > raw_AIMs_list



