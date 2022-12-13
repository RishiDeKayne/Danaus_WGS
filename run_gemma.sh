#!/bin/bash
#this script will run a gwas
#run as follows:
#first run sconda gemma
#./run_gemma.sh plink_output_path_and_prefix phenotypes run_folder

#e.g. ./run_gemma.sh /scratch/rdekayne/gemma/relatedness /scratch/rdekayne/gemma/background_colour_values.csv /scratch/rdekayne/gemma/background

FAM=$1
PHENO=$2
WORK_DIR=$3

cd ${WORK_DIR}

#take our hindwing phenotype file: background_colour_values.csv and add to fam file in 6th column

awk '{print $1,$2,$3,$4,$5}' ${FAM}.fam > relatedness_start.txt

paste relatedness_start.txt ${PHENO} > pre_background.fam

cp ../relatedness/relatedness* .
rm /scratch/rdekayne/gemma/background/relatedness.fam 

awk '{print $1,$2,$3,$4,$5,$7}' ./pre_background.fam > /scratch/rdekayne/gemma/background/relatedness.fam

gemma -bfile relatedness -gk 2 -o gemma_in_rel_background

cp output/* .

cp ../gemma_in* .
rm gemma_in.fam
cp relatedness.fam gemma_in.fam

gemma -bfile gemma_in -k gemma_in_rel_background.sXX.txt -lmm 4 -o background
