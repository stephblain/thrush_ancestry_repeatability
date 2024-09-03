#!/bin/bash

REF="/scratch/user/sblain/get_AIMs/define_AIMs/inland_ref/inland_genome.fasta"
RUN_LOCATION="/scratch/user/sblain/get_AIMs/ref_panel/"
BAM_LOCATION="/scratch/user/sblain/get_AIMs/ref_panel/ref_bams/"
SCAFFOLD_LIST="/scratch/user/sblain/get_AIMs/ref_panel/scaffold_list"
AIMS_PATH="/scratch/user/sblain/get_AIMs/define_AIMs/ancestry_ref_panel_1_m2/raw_AIMs_pos2"
run_tag="ref"


##EDIT sed lines to match number of scaffolds (or which scaffolds, if wanted)
sed -n 1,161p $SCAFFOLD_LIST | while read SCAFFOLD
do
        echo "#!/bin/bash

#SBATCH --job-name="$run_tag"_"$SCAFFOLD"
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=60G
#SBATCH --cpus-per-task=1
#SBATCH --output=error_"$run_tag"_"$SCAFFOLD".%j


module load GCC/11.2.0 BCFtools/1.14

cd $RUN_LOCATION/"$run_tag"_variantCalls

bcftools mpileup -r $SCAFFOLD -o "$SCAFFOLD"."$run_tag".bcf -f $REF *.combo.q30.unbiased.sort.bam

echo "$SCAFFOLD mpileup done"

bcftools call -mO z -o "$SCAFFOLD"."$run_tag".vcf.gz "$SCAFFOLD"."$run_tag".bcf

echo "$SCAFFOLD calling done"

tabix "$SCAFFOLD"."$run_tag".vcf.gz

echo "$SCAFFOLD vcf indexed"

bcftools view -R $AIMS_PATH "$SCAFFOLD"."$run_tag".vcf.gz -o "$SCAFFOLD"."$run_tag".AIMs.vcf.gz

echo "$SCAFFOLD aims selected"

tabix "$SCAFFOLD"."$run_tag".AIMs.vcf.gz

echo "$SCAFFOLD aims vcf indexed"" > vcf_scripts/bams_to_aims_"$SCAFFOLD".sh

#gives permission to write/create a file
chmod 755 vcf_scripts/bams_to_aims_"$SCAFFOLD".sh

#runs script - do a first run with this commented out to make sure
#that writing scripts works
sbatch vcf_scripts/bams_to_aims_"$SCAFFOLD".sh

done