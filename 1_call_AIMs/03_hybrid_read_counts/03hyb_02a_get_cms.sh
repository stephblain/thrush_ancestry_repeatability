#!/bin/bash
      
#SBATCH --job-name=dist_cMs
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --output=e_dist_cMs_1.%j

cd /scratch/user/sblain/get_AIMs/hmm_AIMs/lowCov/get_f50_cMs

module load GCC/9.3.0  OpenMPI/4.0.3 pandas/1.1.2-Python-3.8.2

python ../scripts/03hyb_02b_dist_cMs_byScaf.py

cat current_aims_f50_s* > current_aims_f50.cMs
