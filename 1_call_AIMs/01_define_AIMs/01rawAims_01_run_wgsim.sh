
#!/bin/bash
      
##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=wgsim_coastal   	#Set the job name to "JobExample4"
#SBATCH --time=02:00:00          	#Set the wall clock limit to 1hr and 30min
#SBATCH --ntasks=1               	#Request 1 task
#SBATCH --mem=8G              		#Request 2560MB (2.5GB) per task
#SBATCH --output=error_out.%j  		#Send stdout/err to "Example4Out.[jobID]"
#SBATCH --gres=gpu:1             	#Request 1 GPU per node can be 1 or 2
#SBATCH --partition=gpu          	#Request the GPU partition/queue


###NOTE this amount of time was overkill - it only took 17 mins

#note that wgsim is from samtools
#to load samtools 1.14, need to load gcc 11.2.0

module load GCC/11.2.0
module load SAMtools/1.14


##Example from Quinn
#example with no errors (-e 0), 30X coverage for mapping to domestic cat reference genome (-N 255407230), no mutations (-r 0), no indels (-R 1), and paired-end 150 read (-1 150 -2 150)
#wgsim -e 0 -N 255407230 -r 0 -R 1 -1 150 -2 150 serval_draft_genome_October2018_Murphy.fa serval.read1.fq serval.read2.fq
#gzip serval.read1.fq
#gzip serval.read2.fq

# -N is the number of read pairs, which is coverage*[genome length]/([read length]*2)
# I used 30*1131635156/(150*2)
# the extra *2 in the denominator is because it's read pairs, not reads

wgsim -e 0 -N 113163516 -r 0 -R 1 -1 150 -2 150 /scratch/user/sblain/ref_genomes/coastal_genome.fasta coastal.read1.fastq coastal.read2.fastq