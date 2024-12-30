# Call ancestry informative markers in Swainson's thrushes

Pipeline to get AIMs

Adapted from ancestryInfer

## 01 find divergent markers

* Starting from two fasta sequences for reference genomes, convert one to a pseudoref that is collinear with the other

### 01rawAims_01_run_wgsim.sh

* Simulate fastqs from fasta

### 01rawAims_02_prep_ref.sh & 01rawAims_03_align_sort.sh & 01rawAims_04_call_variants.sh & 01rawAims_05_filter.sh

* Align, call, and filter variants

### 01rawAims_06_pseudoref_to_aims.sh

* Get a consensus sequence / pseudoref
* Identify a list of raw AIMs (all fixed differences between the ref and collinear pseudoref)

## 02 reference panel

### 02refPanel_01_sortBams.sh

* Sort aligned bams

### 02refPanel_02_bams_to_aimsVcf.sh

* Call variants
* Select raw AIMs

### 02refPanel_03_allele_counts.sh

* Estimate allele frequency diffs between species - filter based on cutoff difference

### 02refPanel_04_count_outputs.R & 02refPanel_04a_count_outputs_mtDNA.R & 02refPanel_04b_count_outputs_checks.R

* Estimate allele counts at filtered AIMs

## 03 add the hybrids

### 03hyb_01a_write_variantCall_slurms.py & 03hyb_01b_variantCall.sh

* Process low coverage genomes from bams to read counts

### 03hyb_02a_get_cms.sh & 03hyb_02b_dist_cMs_byScaf.py

* Use a recombination map to estimate the distance in cMs between AIMs

### 03hyb_03_format_for_hmm.sh

* Sort files and combine the parent allele counts, cMs between AIMs, and hybrids

## 04 run the markov model

### 04hmm_01_run_hmm.sh

* Run AncestryHMM

### 04hmm_02a_combine_posteriors.sh

* Get posterior probs of each state

### 04hmm_02b_get_posterior_plinkFormat.r & 04hmm_02b_get_posterior_states.r

* Get estimated ancestry state for sites above cutoff; output in table and plink formats

### 04hmm_03_LDprune.sh

* LD-prune AIMs from plink format
