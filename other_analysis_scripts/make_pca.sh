#!/bin/bash
#first run: sconda PCA
#you may have to run chmod +x make_pca.sh
#usage: ./make_pca.sh my.vcf output_prefix

#this script will produce two types of PCA, the first linkage filtered with the plink parameters --indep-pairwise 50 10 0.1 and the other unfiltered
#important files output are: *_filt_out.eigenval, *_filt_out.eigenvec and *_unfilt_out.eigenval, *_unfilt_out.eigenvec
VCF=$1
PREFIX=$2

#Filtered PCA
plink --vcf ${VCF} --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out ${PREFIX}_filt_out

plink --vcf ${VCF} --double-id --allow-extra-chr --set-missing-var-ids @:# --extract ${PREFIX}_filt_out.prune.in --make-bed --pca --out ${PREFIX}_filt_out

#Unfiltered PCA
plink --vcf ${VCF} --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out ${PREFIX}_unfilt_out

echo "done"
