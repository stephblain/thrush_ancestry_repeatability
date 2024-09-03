######Formatting

module purge
module load GCC/11.2.0 BCFtools/1.14 VCFtools/0.1.16

#vcf of assorted Catharus species
#from Delmore et al.; in review at Molecular Ecology
vcf=genotype_all_fix_ust.vcf.gz

##pull out relevent chunks of the genome with bcftools

chrBed=chr5_block.bed
prefix=catharusChr5Block

bcftools view -R $chrBed -m2 -M2 -v snps $vcf -O z -o $prefix.vcf.gz
tabix $prefix.vcf.gz

chrBed=chr1_block.bed
prefix=catharusChr1Block

bcftools view -R $chrBed -m2 -M2 -v snps $vcf -O z -o $prefix.vcf.gz
tabix $prefix.vcf.gz

#minimum 8 samples to include site
python $SCRATCH/tools/vcf2phylip.py -i $prefix.vcf.gz --output-prefix $prefix -m 8 -o KSW5403

module purge
module load GCC/11.3.0 OpenMPI/4.1.4 IQ-TREE/2.2.2.3


prefix=catharusChr1Block
iqtree2 -s $prefix.min8.phy --prefix $prefix -m MFP #select best model with BIC - GTR+F for Chr1
iqtree2 -s $prefix.min8.phy --prefix $prefix -m GTR+F -B 10000 -redo #fit tree with 10,000 bootstrap replicates


prefix=catharusChr5Block
iqtree2 -s $prefix.min8.phy --prefix $prefix -m MFP
iqtree2 -s $prefix.min8.phy --prefix $prefix -m TVM+F -B 10000 -redo
