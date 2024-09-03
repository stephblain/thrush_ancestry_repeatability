module load GCC/9.3.0 PLINK/1.9b5

################# Autosomal AIMs

prefix="thrushHybrids.f50.20240501" #prefix preceding states, if there is a states.map and states.ped present
statesFile="/scratch/user/sblain/get_AIMs/run_hmm/autosomes_f50/Mar2024_byScaf/thrushHybrids.f50.20240501.states"
probsFile="/scratch/user/sblain/get_AIMs/run_hmm/autosomes_f50/Mar2024_byScaf/thrushHybrids.f50.20240501.probs"

cd /scratch/user/sblain/get_AIMs/run_hmm/autosomes_f50/Mar2024_byScaf



# apply filters - missing data < 0.25 and minor allele frequency > 0.05
plink --file $prefix.states \
	--maf 0.05 --geno 0.25 \
	--allow-extra-chr --indep-pairwise 200 20 .2 \
	--out $prefix.miss25 \
	--set-missing-var-ids @:# \
	--double-id

sed 's/_/\t/g' $prefix.miss25.prune.in > $prefix.miss25.prune.regions #make a tab-delimited file

sed 's/\t/_/' $statesFile > $statesFile.mod #make mod file with chrom_pos
#grep search for chrom_pos then reformat to non-mod file
grep -wFf $prefix.miss25.prune.in $statesFile.mod | \
 sed 's/_/\t/' > $prefix.LD.states

#add in header line with references etc
(head -n 1 $statesFile && cat $prefix.LD.states) \
  > tmp && mv tmp $prefix.LD.states
  
  
sed 's/\t/_/' $probsFile > $probsFile.mod #make mod file with chrom_pos
#grep search for chrom_pos then reformat to non-mod file
grep -wFf $prefix.miss25.prune.in $probsFile.mod | \
 sed 's/_/\t/' > $prefix.LD.probs

#add in header line with references etc
(head -n 1 $probsFile && cat $prefix.LD.probs) \
  > tmp && mv tmp $prefix.LD.probs



################# Z chromosome AIMs - run on males only, then extract same set from females

prefix="thrushHybridsZchromMales.f50.20240429" #prefix preceding states, if there is a states.map and states.ped present
statesFile="/scratch/user/sblain/get_AIMs/run_hmm/ZchromMales_f50/thrushHybridsZchromMales.f50.20240429.states"
probsFile="/scratch/user/sblain/get_AIMs/run_hmm/ZchromMales_f50/thrushHybridsZchromMales.f50.20240429.probs"

cd /scratch/user/sblain/get_AIMs/run_hmm/ZchromMales_f50

# apply filters - missing data < 0.25 and minor allele frequency > 0.05
plink --file $prefix.states \
	--maf 0.05 --geno 0.25 \
	--allow-extra-chr --indep-pairwise 200 20 .2 \
	--out $prefix.miss25 \
	--set-missing-var-ids @:# \
	--double-id

sed 's/_/\t/g' $prefix.miss25.prune.in > $prefix.miss25.prune.regions #make a tab-delimited file

sed 's/\t/_/' $statesFile > $statesFile.mod #make mod file with chrom_pos
#grep search for chrom_pos then reformat to non-mod file
grep -wFf $prefix.miss25.prune.in $statesFile.mod | \
 sed 's/_/\t/' > $prefix.LD.states

#add in header line with references etc
(head -n 1 $statesFile && cat $prefix.LD.states) \
  > tmp && mv tmp $prefix.LD.states
  
  
sed 's/\t/_/' $probsFile > $probsFile.mod #make mod file with chrom_pos
#grep search for chrom_pos then reformat to non-mod file
grep -wFf $prefix.miss25.prune.in $probsFile.mod | \
 sed 's/_/\t/' > $prefix.LD.probs

#add in header line with references etc
(head -n 1 $probsFile && cat $prefix.LD.probs) \
  > tmp && mv tmp $prefix.LD.probs

cp $prefix.miss25.prune.in ../ZchromFemales_f50

#####extract from females dataset
  
prefix="thrushHybridsZchromFemales.f50.20240429" #prefix preceding states, if there is a states.map and states.ped present
statesFile="/scratch/user/sblain/get_AIMs/run_hmm/ZchromFemales_f50/thrushHybridsZchromFemales.f50.20240429.states"
probsFile="/scratch/user/sblain/get_AIMs/run_hmm/ZchromFemales_f50/thrushHybridsZchromFemales.f50.20240429.probs"
pruneIn="thrushHybridsZchromMales.f50.20240429.miss25.prune.in"
  
cd /scratch/user/sblain/get_AIMs/run_hmm/ZchromFemales_f50


sed 's/\t/_/' $statesFile > $statesFile.mod #make mod file with chrom_pos
#grep search for chrom_pos then reformat to non-mod file
grep -wFf $pruneIn $statesFile.mod | \
 sed 's/_/\t/' > $prefix.LD.states

#add in header line with references etc
(head -n 1 $statesFile && cat $prefix.LD.states) \
  > tmp && mv tmp $prefix.LD.states
  
  
sed 's/\t/_/' $probsFile > $probsFile.mod #make mod file with chrom_pos
#grep search for chrom_pos then reformat to non-mod file
grep -wFf $pruneIn $probsFile.mod | \
 sed 's/_/\t/' > $prefix.LD.probs

#add in header line with references etc
(head -n 1 $probsFile && cat $prefix.LD.probs) \
  > tmp && mv tmp $prefix.LD.probs