import pandas as pd

#Bams to vcf for high coverage sequences from the thrush family
#Hermit, Bicknell's, Quebec inland Swainson's, BC inland Swainson's, BC coastal Swainson's

samples_df = pd.read_table('sample_list_dstats.tsv').set_index("sample", drop=False)
SAMPLES= list(samples_df['sample'])

scaffolds_df = pd.read_table('scaffold_list.tsv').set_index("scaffold", drop=False)
SCAFFOLDS= list(scaffolds_df['scaffold'])


REF="/scratch/user/sblain/get_AIMs/define_AIMs/inland_ref/inland_genome.fasta"

#module load GCC/11.3.0 SAMtools/1.16.1 OpenMPI/4.1.4 snakemake/7.22.0
#run as:  nohup snakemake --profile . &
#before running, make config.yaml file in directory

###Before running, format bam files names to be: SAMPLENAME.combo.q30.unbiased.bam with an associated .bai
## put bams in ./dstat_bams/ and make a directory called ./dstat_calls/

rule all:
	input:
		"thrush.d_stats.filtered.vcf.gz.tbi"

rule sort_bam:
	input:
		bam=expand("dstat_bams/{{sample}}.combo.q30.unbiased.bam", sample=SAMPLES),
		bai=expand("dstat_bams/{{sample}}.combo.q30.unbiased.bam.bai", sample=SAMPLES)
	output:
		bam=expand("dstat_bams/{{sample}}.combo.sort.bam", sample=SAMPLES)		
	shell:
		"samtools sort {input.bam} -o {output.bam}"
		

rule index_sorted_bams:
	input:
		bam=expand("dstat_bams/{{sample}}.combo.sort.bam", sample=SAMPLES)
	output:
		bai=expand("dstat_bams/{{sample}}.combo.sort.bam.bai", sample=SAMPLES)
	shell:
		"samtools index {input.bam}"



rule bcftools_call:
	input:
		fa=REF,
		bam=expand("dstat_bams/{sample}.combo.sort.bam", sample=SAMPLES),
		bai=expand("dstat_bams/{sample}.combo.sort.bam.bai", sample=SAMPLES),
	params:
		scaf=expand("{{scaffold}}", scaffold=SCAFFOLDS)
	threads: 4
	output:
		bcf=expand("dstat_calls/thrush.refOut.{{scaffold}}.bcf",scaffold=SCAFFOLDS)
	shell:
		"bcftools mpileup -r {params.scaf} -f {input.fa} {input.bam} | "
		"bcftools call -mO b -o {output.bcf}"

#concatenate bcf files
#keep only biallelic SNPs
#filter out sites with > 25% missing data or a minor allele frequency < 0.05
rule bcftools_concat:
	input:
		bcf=expand("dstat_calls/thrush.refOut.{scaffold}.bcf",scaffold=SCAFFOLDS)
	output:
		vcfCat="thrush.d_stats.filtered.vcf.gz"
	shell:
		"bcftools concat {input.bcf} | "
		"bcftools view -m2 -M2 -v snps | "
		"bcftools filter -e 'F_MISSING > 0.25 || MAF <= 0.05' -O z -o {output.vcfCat}"

rule index_vcf:
	input:
		vcf="thrush.d_stats.filtered.vcf.gz"
	output:
		tbi="thrush.d_stats.filtered.vcf.gz.tbi"
	shell:
		"tabix {input.vcf}"

