#!/bin/bash

#SBATCH --job-name=align_sort   	#Set the job name
#SBATCH --time=12:00:00          	#Set the wall clock limit
#SBATCH --ntasks=3                      #Request 3 tasks
#SBATCH --mem=8G              		#Request 8GB per task
#SBATCH --output=error_03_align_sort.%j	#Send stdout/err to doc
#SBATCH --gres=gpu:1             	#Request 1 GPU per node can be 1 or 2
#SBATCH --partition=gpu          	#Request the GPU partition/queue

cd /scratch/user/sblain/get_AIMs/define_AIMs

module load GCCcore/11.2.0
module load BWA/0.7.17
module load GCC/11.2.0
module load SAMtools/1.14


#align with bwa-mem
bwa mem -M /scratch/user/sblain/ref_genomes/inland_genome.fasta \
  -t 3 \
  coastal.read1.fastq coastal.read2.fastq > coastal.samOut.sam

#sort, filter, and index with samtools

#fills in info on corresponding paired end reads
samtools fixmate -O bam coastal.samOut.sam coastal.bamOut.bam
#sort reads and index
samtools sort coastal.bamOut.bam -o coastal.bamOut.sorted.bam
samtools index coastal.bamOut.sorted.bam
#filters out qual less than 30, output in bam
samtools view -b -q 30 coastal.bamOut.sorted.bam > coastal.bamOut.sorted.unique.bam 
samtools index coastal.bamOut.sorted.unique.bam


