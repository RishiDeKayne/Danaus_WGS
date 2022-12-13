#!/bin/bash
#this script will run a gwas
#run as follows:
#first run sconda gemma
#./run_gemma.sh relatedness_output_path_and_prefix gemma_in_path_and_prefix phenotypes run_folder name

#e.g. ./run_gemma.sh /scratch/rdekayne/gemma/relatedness /scratch/rdekayne/gemma/gemma_in /scratch/rdekayne/gemma/background_colour_values.csv /scratch/rdekayne/gemma/background rel_background

REL=$1
FAM=$2
PHENO=$3
WORK_DIR=$4
OUT_NAME=$5

cd ${WORK_DIR}

#take our hindwing phenotype file: background_colour_values.csv and add to fam file in 6th column

awk '{print $1,$2,$3,$4,$5}' ${REL}.fam > relatedness_start.txt

paste relatedness_start.txt ${PHENO} > pre_background.fam

cp ${REL}* .
rm -f ${WORK_DIR}/relatedness.fam 

awk '{print $1,$2,$3,$4,$5,$7}' ./pre_background.fam > ${WORK_DIR}/relatedness.fam

gemma -bfile relatedness -gk 2 -o ${OUT_NAME}

cp output/* .

cp ${FAM}* .
rm gemma_in.fam
cp relatedness.fam gemma_in.fam

gemma -bfile gemma_in -k ${OUT_NAME}.sXX.txt -lmm 4 -o ${OUT_NAME}_background_full
