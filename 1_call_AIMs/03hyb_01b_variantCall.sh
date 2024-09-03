#!/bin/bash

#run as:
#bash 01_bams_to_counts.sh REF BAM_LOCATION IND_LIST RENAME_SCAFFOLDS RUN_LOCATION run_tag n1 n2 AIMS_LIST A_INFER RECOMB_RATE
#see variable descriptions in "01_write_variantCall_slurms.py"


#input data: bams that have been aligned to both subspecies genomes
#this script sorts bams, uses bcftools to call variants (mpileup + call), rewrites scaffold names in vcf to replace underscores with dashes, then deletes intermediate files

#NOTE: this script relies on bam files having the name format "indivdual ID".combo.q30.unbiased.bam (ex. AH28H09.combo.q30.unbiased.bam)



#take variables from command line input
REF=$1
BAM_LOCATION=$2
IND_LIST=$3
RENAME_SCAFFOLDS=$4
RUN_LOCATION=$5
run_tag=$6
n1=$7
n2=$8
AIMS_LIST=$9
A_INFER=${10}
RECOMB_RATE=${11}

cd $RUN_LOCATION

if [ ! -d "$run_tag"_variantCalls ]
then
	mkdir "$run_tag"_variantCalls
fi

cd "$run_tag"_variantCalls



sed -n "$n1","$n2"p $IND_LIST | while read IND

	if [ -z $IND ] #avoid infinite loop if variable is empty
	then
		echo "individual variable empty"
		break
	fi
	
	do
		if [[ -s "$IND".noUnderscores.vcf || -s "$IND".noUnderscores.vcf_counts.hmm ]]
		#only run through variant calling etc if vcf or hmm file is not empty
		then
			echo "$IND vcf or hmm file exists"
		else
			#first sort bams
			
			SORT_BAM="$IND".combo.q30.unbiased.sort.bam
			
			if [[ -s "$SORT_BAM" ]] #check if sort bam file exists
			then
				echo "$SORT_BAM exists. Not re-sorting bam."
			else

				BAM="$BAM_LOCATION""$IND".combo.q30.unbiased.bam
				samtools sort $BAM -o $SORT_BAM
				echo "$IND bam sorted"
			fi	
			
			#then variant calling
			
			MPILEUP="$IND".bcf
			VCF="$IND".vcf.gz
			
			if [[ -s "$VCF" ]]
			then
				echo "$VCF exists - $MPILEUP finished creation. Not creating bcf again."
			else #variant calling step one
				bcftools mpileup -o $MPILEUP -f $REF "$IND".combo.q30.unbiased.sort.bam
				echo "$IND mpileup bcf created"
			fi
			
				
			if [[ -s "$VCF".tbi ]]
			then
				echo "$VCF and tbi exist. Not re-calling snps"
			else #variant calling step two
				bcftools call -mO z -o $VCF $MPILEUP
				tabix $VCF
				echo "$IND variants called"
			fi
		fi
		
		#then convert to no underscores format
		#This block of code could be more efficient - maybe with a while loop
				
		if [[ -s "$IND".noUnderscores.vcf ]]
		#check if vcf exists and is not empty
		then
			
			#check if last line of vcf includes "mtDNA"
			mtLine=$(tail -1 "$IND".noUnderscores.vcf | cut -c1-5)
			scaf9Line=$(tail -1 "$IND".noUnderscores.vcf | cut -c1-11)
			
			if [[ $mtLine == "mtDNA" || $scaf9Line == "scaffold-9-" ]]
			then
				echo "$IND vcf already converted to no underscores format - end in $scaf9Line"
			else #rename scaffolds in vcf to have dashes instead of underscores - necessary for ancestryInfer	
				bcftools annotate --rename-chrs $RENAME_SCAFFOLDS "$IND".vcf.gz > "$IND".noUnderscores.vcf
				echo "$IND scaffolds renamed"
			fi
		elif [[ -s "$IND".noUnderscores.vcf_counts.hmm ]]
		#if vcf is empty but the hmm file exists
		then
			echo "$IND hmm already exists."
		else #rename scaffolds in vcf to have dashes instead of underscores - necessary for ancestryInfer
			bcftools annotate --rename-chrs $RENAME_SCAFFOLDS "$IND".vcf.gz > "$IND".noUnderscores.vcf
			echo "$IND scaffolds renamed"
		fi
		
		#then delete intermediate files
		
		if [[ -s "$IND".noUnderscores.vcf && ! -s "$IND".noUnderscores.vcf_counts.hmm ]]
		#only if vcf not empty and hmm (end output) does not exist
		then
		
			#check if last line of vcf includes "mtDNA" or "scaffold-9-" before deleting intermediate files
			mtLine1=$(tail -1 "$IND".noUnderscores.vcf | cut -c1-5)
			scaf9Line1=$(tail -1 "$IND".noUnderscores.vcf | cut -c1-11)
			if [[ $mtLine1 == "mtDNA" || $scaf9Line1 == "scaffold-9-" ]]
			then
				rm "$IND".combo.q30.unbiased.sort.bam "$IND".bcf "$IND".vcf.gz "$IND".vcf.gz.tbi
				echo "$IND intermediate files deleted after variant calling - last line was $scaf9Line"
			else
				echo "$IND intermediate files not deleted after variant calling - check variant calls"
			fi
		else
			echo "$IND vcf doesn't exist and or hmm does exist"
		fi
		
		#check that vcf exists and bam file was deleted
		if [[ -s "$IND".noUnderscores.vcf && ! -a "$IND".combo.q30.unbiased.sort.bam ]]
		then
			
			echo "$IND vcf exists. Start read counting."
			
			VCF="$IND".noUnderscores.vcf
			perl "$A_INFER"/vcf_to_counts_non-colinear.pl $VCF $AIMS_LIST $A_INFER
			#expected output: _counts, .aims, and .mod files
			

			#then format read counts for HMM
			#if aims file was created properly
			
			if [[ -s "$VCF".aims ]]
			then
		
				#define variables
				COUNTS1="$VCF"_counts
				HMMSITES1="$COUNTS1".hmm
				COUNTS_PRESENT="$COUNTS1"_present
				COUNTS_ABSENT="$COUNTS1"_absent
				COUNTS_COMBINED="$COUNTS1"_combined
				COUNTS_COMBINED_SORTED="$COUNTS1"_combined_sorted
				COUNTS_BED="$COUNTS1".bed
							
				#reformat into bed file
				cat $COUNTS1 | perl -p -e 's/_/\t/g' | awk -v OFS="\t" \$1=\$1\"\\t\"\$2\ > $COUNTS_BED

				AIMS_BED="$AIMS_LIST".mod.bed

				#find overlap (and then difference) between AIMs bed file and individual bed file
				bedtools intersect -a $AIMS_BED -b $COUNTS_BED -wb -f 1 > $COUNTS_PRESENT
				bedtools intersect -a $AIMS_BED -b $COUNTS_BED -wb -v -f 1 | awk -v OFS="\t" \$6\=\$1\"\\t\"\$2\"\\t\"\$2\"\\t\"\$4\"\\t\"\$5\"\\t0\\t0\" > $COUNTS_ABSENT

				#make one big happy bed file
				cat $COUNTS_PRESENT $COUNTS_ABSENT > $COUNTS_COMBINED

				#run ancestry infer perl scripts to reformat into hmm input files and hope the output is correct
				perl $A_INFER/merge_files_using_two_columns_sharing_values_stdout.pl $AIMS_BED 0 1 $COUNTS_COMBINED 0 1 | cut -f 1-5 --complement > $COUNTS_COMBINED_SORTED
				perl $A_INFER/reformat_counts_v2.pl $COUNTS_COMBINED_SORTED $RECOMB_RATE > $HMMSITES1
				
				echo "$IND reads counted."
			else
				echo "$IND aims file empty or missing."
			fi
			#divide into individual and parent counts by splitting columns
			
			# NEXT: only make parental files once and use updated read counts for that
			
			# HMM="$IND".hmm.combined
			# PAR_CURR="$HMM".parental.format
			# IND_CURR="$HMM".pass.formatted

			# cut -f 1-6 $HMMSITES1 > $PAR_CURR #no need to make a parental format file for every individual
			# cut -f 7-8 $HMMSITES1 > $IND_CURR
			
			# make list of read count files
			# echo "$IND_CURR" >> ../HMM.hybrid.files.list_allchrs
			# echo "$PAR_CURR" >> ../HMM.parental.files.list_allchrs
		fi
		
		if [[ -s $COUNTS_PRESENT && -s $HMMSITES1 ]] #if counts are actually there and pipeline made it to last step
		then
		
			#check if last line of vcf includes "mtDNA" before deleting intermediate files
			mtLine2=$(tail -1 $HMMSITES1 | cut -c1-5)
			
			if [[ $mtLine2 == "mtDNA" ]]
			then
				rm "$IND".noUnderscores.vcf_counts "$IND".noUnderscores.vcf_counts_combined_sorted "$IND".noUnderscores.vcf_counts_absent
				rm "$IND".noUnderscores.vcf  "$IND".noUnderscores.vcf_counts.bed "$IND".noUnderscores.vcf.mod
				rm "$IND".noUnderscores.vcf_counts_present "$IND".noUnderscores.vcf.aims "$IND".noUnderscores.vcf_counts_combined
				echo "$IND intermediate counts files deleted after read counting"
			else
				echo "$IND intermediate files not deleted after read counting - check count estimation"
			fi
		fi

	done
