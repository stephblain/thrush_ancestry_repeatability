# Call ancestry informative markers in Swainson's thrushes

Pipeline to get AIMs

Adapted from ancestryInfer

## 01 find divergent markers

Starting from two fasta sequences for reference genomes, convert one to a pseudoref that is collinear with the other

Simulate fastqs from fasta

Align and filter variants

Get a consensus sequence / pseudoref

Identify a list of raw AIMs (all fixed differences between the ref and collinear pseudoref)

## 02 reference panel

Call variants

Select raw AIMs

Estimate allele frequency diffs between species in ref panel 1 - filter based on cutoff difference

Estimate allele counts at filtered AIMs in ref panel 2

## 03 add the hybrids

Process low coverage genomes from bams to read counts

Use a recombination map to estimate the distance in cMs between AIMs

Sort files and combine the parent allele counts, cMs between AIMs, and hybrids

## 04 run the markov model

Run AncestryHMM

Get posterior probs of each state

Get estimated ancestry state for sites above cutoff; output in table and plink formats
