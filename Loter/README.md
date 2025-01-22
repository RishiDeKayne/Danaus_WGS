
### Scripts and data for making and plotting Loter results.

The script `loter_wrapper.py` was used to run loter as follows:

```bash
python loter_wrapper.py -i contig15.1.double.phased.minVar2HET60BI.haplo_count.tsv.gz -o contig15.1.double.phased.minVar2HET60BI.loter_K4_160noext.tsv.gz -g A SM21TS03,SM21TS05,SM21TS16,SM21TS33,SM21TS34,SM21TS35,SM21TS36 -g B SM15W61,SM15W66,SM15W69,SM15W72,SM15W74,SM18W01 -g C SM16N01,SM16N04,SM16N05,SM16N06,SM16N20,SM16N37 -g D SM21SP06,SM21SP09 --haploidify_IDs --min_score 160 --threads 10
```

Note that the option --min_score 160 means that the script only assigns ancestry when the score is 160/160. This was used for the chromosome-wide view shown in Supplementary Text S1 figure J. A lower threshold of 140/160 was used for the close-up ancestry pattern shown in Figure 4B and 4D.

Then [`phasepaint`](https://github.com/simonhmartin/phasepaint) was used to correct obvious phasing errors

```bash
python phasepaint.py -i contig15.1.double.phased.minVar2HET60BI.loter_K4_160noext.tsv.gz -o contig15.1.double.phased.minVar2HET60BI.loter_K4_160noext.phasepaint_I20.tsv.gz --threads 20 --max_iterations 20
```
