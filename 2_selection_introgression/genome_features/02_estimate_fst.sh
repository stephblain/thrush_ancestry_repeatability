# get average fst in sliding 100kb windows

VCF=all_scaffolds.filtered.ref.vcf.gz
pop1=inds_ref_inland
pop2=inds_ref_coastal

cd /scratch/user/sblain/get_AIMs/ref_panel

module load  GCC/11.2.0 BCFtools/1.14 VCFtools/0.1.16 #module changeover for vcftools

# estimate fst
vcftools --gzvcf $VCF --weir-fst-pop $pop1 --weir-fst-pop $pop2 --fst-window-size 100000 --fst-window-step 100000 --out thrushRef

