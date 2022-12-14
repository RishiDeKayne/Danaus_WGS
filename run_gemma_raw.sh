#!/bin/bash
#this script will run a gwas
#run as follows:
#first run sconda gemma
#./run_gemma.sh gemma_in_path_and_prefix phenotypes run_folder name

#e.g. ./run_gemma_raw.sh /scratch/rdekayne/gemma/gemma_in /scratch/rdekayne/gemma/background_colour_values.csv /scratch/rdekayne/gemma/background_raw raw_background

FAM=$1
PHENO=$2
WORK_DIR=$3
OUT_NAME=$4

cd ${WORK_DIR}

awk '{print $1,$2,$3,$4,$5}' ${FAM}.fam > in_start.txt
paste in_start.txt ${PHENO} > pre_background.fam

cp ${REL}* .
rm gemma_in.fam

awk '{print $1,$2,$3,$4,$5,$7}' ./pre_background.fam > gemma_in.fam

gemma -bfile gemma_in -gk 2 -o ${OUT_NAME}

cp output/* .

gemma -bfile gemma_in -k ${OUT_NAME}.sXX.txt -lmm 4 -o ${OUT_NAME}

echo ${OUT_NAME} > summary.${OUT_NAME}.txt
grep "analyzed" output/${OUT_NAME}.log.txt >> summary.${OUT_NAME}.txt
echo 'highest snp' >> summary.${OUT_NAME}.txt
min=`awk 'BEGIN{a=1000}{if ($15<0+a) a=$15} END{print a}' output/${OUT_NAME}_background_full.assoc.txt` ; echo $min >> summary.${OUT_NAME}.txt
