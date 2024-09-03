#!/usr/bin/env python3

#writes a set of scripts that can be run in parallel for variant calling individuals with bcftools
#each script then calls 01_variantCall.sh for actual variant calling

#set vars
run_tag="lowCov" #this will be attached to various files and folders
REF="/scratch/user/sblain/get_AIMs/define_AIMs/inland_ref/inland_genome.fasta" #location of ref genome
BAM_LOCATION="/scratch/user/sblain/get_AIMs/hmm_AIMs/lowCov_bams/" #all bams must be in one folder - give the file path to it
IND_LIST="/scratch/user/sblain/get_AIMs/hmm_AIMs/inds_lowCov_2023" #a list of individuals, one per line
RENAME_SCAFFOLDS="/scratch/user/sblain/get_AIMs/hmm_AIMs/test_run_1/rename_scaffolds.txt" #a file for renaming scaffolds from names with underscores to names with dashes. One scaffold per line with old name then new name, tab-delimited
RUN_LOCATION="/scratch/user/sblain/get_AIMs/hmm_AIMs/lowCov/" #where you want output folder to be created
LOAD_MODULES="module load GCC/11.2.0 BCFtools/1.14" #modules to load - need bcftools
SCRIPT_LOCATION="/scratch/user/sblain/get_AIMs/hmm_AIMs/lowCov/scripts/" #where 01_variantCall.sh lives
INDS_PER=30 #inds per run
AIMS_LIST="/scratch/user/sblain/get_AIMs/hmm_AIMs/test_run_1/current_aims_file"
A_INFER="/scratch/user/sblain/ancestryinfer"
RECOMB_RATE=0.00014
SLURM_SPECS="#SBATCH --time=24:00:00\n#SBATCH --ntasks=4\n#SBATCH --mem=200G"

#Run from the directory where you want the set of scripts to be created

with open(IND_LIST, 'r') as indFile:
    N_LINES = len(indFile.readlines()) #read the number of lines from the individuals file
#N_LINES=564 #or write in the number of individuals directly

N_FILES=int(N_LINES/INDS_PER)+1

for i in range(0,N_FILES):
    
    #make new file and write file header
    f = open("01_bams_to_counts_"+run_tag+"_"+str(i)+".slurm","w")
    
    #edit slurm submission variables here as needed
    f.write("#!/bin/bash\n#SBATCH --job-name=call_variants_"+str(i)+"_"+run_tag+"\n"+SLURM_SPECS+"\n#SBATCH --output=error_01_call_variants_"+str(i)+"_"+run_tag+".%j\n\n")
    f.close()
    
    #define 
    n1=i*INDS_PER+1 #starting individual index
    n2=n1+INDS_PER-1 #ending individual index
    
    #append to existing file
    f = open("01_bams_to_counts_"+run_tag+"_"+str(i)+".slurm","a")
    f.write(LOAD_MODULES+"\n\n")
    f.write("cd "+SCRIPT_LOCATION+"\n\n")
    f.write("bash 01_bams_to_counts.sh "+REF+" "+BAM_LOCATION+" "+IND_LIST+" "+RENAME_SCAFFOLDS+" "+RUN_LOCATION+" "+run_tag+" "+str(n1)+" "+str(n2)+" "+AIMS_LIST+" "+A_INFER+" "+str(RECOMB_RATE)+"\n\n")
    f.close()

#

