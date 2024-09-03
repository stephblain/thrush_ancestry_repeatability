#!/bin/bash

ref_path=/scratch/user/sblain/get_AIMs/ref_panel
RAiSD_path=/scratch/user/sblain/tools/RAiSD/raisd-master/bin/release
run_path=/scratch/user/sblain/thrush_hybrids/positiveSelection

cd $run_path

module purge
module load GCC/13.2.0 GSL/2.7

$RAiSD_path/RAiSD -n coastal -I $ref_path/all_scaffolds.filtered.coastalRef.vcf

echo "end coastal run"

$RAiSD_path/RAiSD -n inland -I $ref_path/all_scaffolds.filtered.inlandRef.vcf

echo "end inland run"