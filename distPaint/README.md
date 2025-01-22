
### Data and R script for ancestry painting plots in Figure 2C and Supplementary File S1 Figures I and K.

Ancestry painting was performed using `distPaint.py` available at https://github.com/simonhmartin/genomics_general.
This script a assigns ancestry in windows along the genome based on dxy between each haplotype and a defined set of reference haplotypes.

The command used was:

```bash
python genomics_general/distPaint.py -w 200 --windType sites -g contig15.1.double.phased.haplo.gz -o contig15.1.double.phased.distPaint.w200p01.tsv.gz -p klugii SM15W61_A,SM15W66_A,SM15W69_A,SM15W72_A,SM15W74_A,SM18W01_A,SM15W61_B,SM15W66_B,SM15W69_B,SM15W72_B,SM15W74_B,SM18W01_B -p orientis SM21TS03_A,SM21TS05_A,SM21TS16_A,SM21TS33_A,SM21TS34_A,SM21TS35_A,SM21TS36_A,SM21TS03_B,SM21TS05_B,SM21TS16_B,SM21TS33_B,SM21TS34_B,SM21TS35_B,SM21TS36_B  -p chrysippus SM16N01_A,SM16N04_A,SM16N05_A,SM16N06_A,SM16N20_A,SM16N37_A,SM16N01_B,SM16N04_B,SM16N05_B,SM16N06_B,SM16N20_B,SM16N37_B -p karamu SM21SP06_A,SM21SP06_B,SM21SP09_A,SM21SP09_B --p_threshold 0.01 --writeFailedWindows -T 20
```

Then [`phasepaint`](https://github.com/simonhmartin/phasepaint) was used to correct obvious phasing errors

```bash
python ~/Research/dev/phasepaint/phasepaint.py -i contig15.1.double.phased.distPaint.w200p01.tsv.gz -o contig15.1.double.phased.distPaint.w200p01.phasepaint_I500.tsv.gz --threads 39 --max_iterations 500 --ignore_first_n_columns 5
```
