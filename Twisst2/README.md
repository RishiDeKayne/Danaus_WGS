
### Data and scripts for Twisst2 analyses.

[`twisst2`](https://github.com/simonhmartin/twisst2) version 0.0.3 was run with the following commands:

#prepare input data
```bash
sticcs prep -i DC174OG.contig15.1.realigned.DP7GQ30het75min87minVar2.vcf.gz --outgroup RG28754 | bgzip > DC174OG.contig15.1.realigned.DP7GQ30het75min87minVar2.DC.vcf.gz
```

#Run twisst2 with three different sets of groups

```bash
twisst2 -i DC174OG.contig15.1.realigned.DP7GQ30het75min87minVar2.DC.vcf.gz --ploidy 2 --max_iterations 100 -o DC174OG.contig15.1.realigned.DP7GQ30het75min87minVar2.Twisst2_K3 --group_names klugii chrysippus orientis --groups  SM15W61,SM15W66,SM15W69,SM15W72,SM15W74,SM18W01 SM16N01,SM16N04,SM16N05,SM16N06,SM16N20,SM16N37 SM21TS03,SM21TS05,SM21TS16,SM21TS33,SM21TS34,SM21TS35,SM21TS36

twisst2 -i DC174OG.contig15.1.realigned.DP7GQ30het75min87minVar2.DC.vcf.gz --ploidy 2 --max_iterations 100 -o DC174OG.contig15.1.realigned.DP7GQ30het75min87minVar2.Twisst2_kar_klu_ori --group_names karamu klugii orientis --groups SM21SP06,SM21SP09 SM15W61,SM15W66,SM15W69,SM15W72,SM15W74,SM18W01 SM21TS03,SM21TS05,SM21TS16,SM21TS33,SM21TS34,SM21TS35,SM21TS36

twisst2 -i DC174OG.contig15.1.realigned.DP7GQ30het75min87minVar2.DC.vcf.gz --ploidy 2 --max_iterations 100 -o DC174OG.contig15.1.realigned.DP7GQ30het75min87minVar2.Twisst2_chryD_chry_ori --group_names chrysippus_dark chrysippus orientis --groups RV12N317,SM19SY01 SM16N01,SM16N04,SM16N05,SM16N06,SM16N20,SM16N37 SM21TS03,SM21TS05,SM21TS16,SM21TS33,SM21TS34,SM21TS35,SM21TS36
```
