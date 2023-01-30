cd /scratch/rdekayne/whatshap_DC174

cat x04_scaffs.txt | while read scaf
do
	touch whatshap_run_$scaf.txt
	cat indiv.list | while read indiv
do
	echo "${indiv}"
	INDIV_BAM=`cat "${indiv}".bam.location.txt`
	bcftools view -s "${indiv}" -r "${scaf}" -o ./"${indiv}"."${scaf}".vcf /data/martin/genomics/analyses/Danaus_popgen/DC174/phased_vcfs/DC174_allchroms_output_filt_mindepth7_minqual30_mac4_AN313_nosmallcontigs.b.vcf.gz 
	bgzip -c ./"${indiv}"."${scaf}".vcf > "${indiv}"."${scaf}".vcf.gz && rm -f "${indiv}"."${scaf}".vcf
	tabix ./"${indiv}"."${scaf}".vcf.gz
	samtools view -b /data/martin/genomics/analyses/Danaus_mapping/DC174/realigned/realigned_bams/"${INDIV_BAM}" "${scaf}" > "${INDIV_BAM}"."${scaf}".bam
	samtools sort "${INDIV_BAM}"."${scaf}".bam > "${INDIV_BAM}"."${scaf}".sorted.bam && samtools index "${INDIV_BAM}"."${scaf}".sorted.bam
	rm -f "${INDIV_BAM}"."${scaf}".bam
	echo 'whatshap phase -o ./phased/'"${indiv}"'.'"${scaf}"'.phased.vcf --reference=/data/martin/genomics/analyses/Danaus_genome/Dchry2/Dchry2.2.fa ./'"${indiv}"'.'"${scaf}"'.vcf.gz '"${INDIV_BAM}"'.'"${scaf}"'.sorted.bam --ignore-read-groups --chromosome '"${scaf}"' --sample='"${indiv}"' && touch ./done/phased.'"${indiv}"'.'"${scaf}"'.done' >> whatshap_run_$scaf.txt
done
done
