# Selection & introgression analysis

## genomic_clines

### 01_run_bgchm_thrushes_withZ.R

* Fit genomic clines using bgc-hm

### 02_process_bgchm_output.Rmd

* Estimate repeatability of cline centre & gradient outliers
* Look for overlap with migratory loci
* Fit models of relationship between ancestry transitions & genome features

## ancestry_transitions

### 01_extract_adult_hybrids.R

* Pre-process data

### 02_find_ancestry_tracts.R

* Find the locations of ancestry transitions based on AIMs genotypes states file

### 03_process_ancestry_transitions.Rmd

* Estimate repeatability of ancestry transitions, look for overlap with migratory loci
* Fit models of relationships between ancestry transitions and genome features

### 04_run_iqtree.sh

* Estimate phylogeny of regions with few transitions

## fd_introgression

### 01_Snakefile.py & 02_format_Dsuite_input.sh

* Put together vcf with relevant SNPs for looking for introgression between coastal & inland BC Swainson's thrushes

### 03_run_Dsuite.sh

* Estimate introgression as D & fd

### 04_process_fd_windows.Rmd

* Fit models of relationships between introgresion and genome features

## genome_features

### 01_estimate_gene_density.R & 02_estimate_fst.sh & 03_estimate_positive_selection.sh & 04_process_positive_selection.ipynb

* Estimate relevant features (Fst, positive selection, gene density) in windows

### 05_combine_genome_features.Rmd

* Pre-process data to combine windowed estimates of genomic features into one dataframe

### 06_plot_genome_features.Rmd

* Make plots of relationships of different variables to genome features

## phenotypes

### phenotypic_repeatability.R

* Estimate relationships between migratory, morphological, colour traits and genomic ancestry
* Calculate repeatability among populations for each trait category & genomic cline gradients