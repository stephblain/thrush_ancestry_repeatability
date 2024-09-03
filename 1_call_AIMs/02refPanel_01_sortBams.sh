#!/bin/bash

#SBATCH --job-name=sortBams
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --output=error_01_sortBams.%j
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu

REF="/scratch/user/sblain/get_AIMs/define_AIMs/inland_ref/inland_genome.fasta"
BAM_LOCATION="/scratch/user/sblain/get_AIMs/ref_panel/ref_bams/"

IND_LIST="/scratch/user/sblain/get_AIMs/ref_panel/ind_list_ref"

module load GCC/11.3.0
module load SAMtools/1.16.1

#cd $BAM_LOCATION
cd "/scratch/user/sblain/get_AIMs/ref_panel/ref_variantCalls"

#first sort bams
while read IND

        if [ -z $IND ]
        #avoid infinite loop if variable is empty
        then
                echo "individual variable empty"
                break
        fi

        do
                BAM="$BAM_LOCATION""$IND".combo.q30.unbiased.bam
                SORT_BAM="$IND".combo.q30.unbiased.sort.bam

                if [[ -f "$SORT_BAM" ]]
                then
                        echo "$SORT_BAM exists. Not re-sorting bam."
                else
                        samtools sort $BAM -o $SORT_BAM
                fi
                samtools index $SORT_BAM

        done < $IND_LIST