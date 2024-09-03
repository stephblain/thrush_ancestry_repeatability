#!/bin/bash
      
##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=prep_inland   	#Set the job name to "JobExample4"
#SBATCH --time=02:00:00          	#Set the wall clock limit to 2 hrs
#SBATCH --ntasks=1               	#Request 1 task
#SBATCH --mem=8G              		#Request 2560MB (2.5GB) per task
#SBATCH --output=error_out.%j  		#Send stdout/err to "Example4Out.[jobID]"
#SBATCH --gres=gpu:1             	#Request 1 GPU per node can be 1 or 2
#SBATCH --partition=gpu          	#Request the GPU partition/queue



####prep the main reference genome, which the pseudo-reference fastq's will be aligned to

cd /scratch/user/sblain/ref_genomes/ #move directories

#load modules

module load picard/2.25.1-Java-11
module load GCCcore/11.2.0
module load BWA/0.7.17
module load GCC/11.2.0
module load SAMtools/1.14

#index genome
bwa index inland_genome.fasta

#make picard dictionary
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary \
  R=inland_genome.fasta O=inland_genome.dict

#index genome a different way
samtools faidx inland_genome.fasta