#!/bin/bash
      
##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=align_sort   	#Set the job name to "JobExample4"
#SBATCH --time=10:00:00          	#Set the wall clock limit to 10 hrs
#SBATCH --ntasks=1               	#Request 1 task
#SBATCH --mem=12G              		#Request 12GB per task
#SBATCH --output=error_04_call_vars.%j	#Send stdout/err to "Example4Out.[jobID]"
#SBATCH --gres=gpu:1             	#Request 1 GPU per node can be 1 or 2
#SBATCH --partition=gpu          	#Request the GPU partition/queue

cd /scratch/user/sblain/get_AIMs/define_AIMs

module load GCC/11.2.0
module load BCFtools/1.14


#generates genotype likelihoods at each position
bcftools mpileup -o coastal.bamOut.sorted.unique.pileup \
  -f /scratch/user/sblain/ref_genomes/inland_genome.fasta \
  coastal.bamOut.sorted.unique.bam

#actually call variants
bcftools call -mO z -o coastal.bamOut.sorted.unique.vcf.gz coastal.bamOut.sorted.unique.pileup
#index called variants
bcftools index coastal.bamOut.sorted.unique.vcf.gz