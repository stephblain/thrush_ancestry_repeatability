#!/bin/bash

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=run_hmm_lowCov
#SBATCH --time=48:00:00
#SBATCH --ntasks=6
#SBATCH --mem=360G
#SBATCH --output=error_run_hmm_lowCov.%j


source activate AncestryEnv
source activate py2 #this one after so bamutils runs
module load GCC/10.3.0 OpenMPI/4.1.1 Armadillo/10.7.5 GSL/2.7

#2 parent populations, with frequency in hybrid samples of 0.78 pop 1 (inland) and 0.22 pop 2 (coastal)
#the minus sign tells the program to estimate things
#start with 3000 generations time estimate for hybridization
ancestry_hmm -a 2 0.78 0.22 -p 0 -10000 -0.78 -p 1 -3000 -0.22 -s lowCov_samples -i f50.hmm.counts.filtered