
### Commands for aligning illumina reads to the Dchry2.2 refernece genome.

The python script `WGS_mapping_bwa_mem.py` is a wrapper for `bwa-mem`, `samtools` and `picard tools` for aliging, sorting and removing duplicate reads.

All read locations are specified in the file `DC174_reads.tsv`

The python wrapper script is run with the following command:

```bash
python WGS_mapping_bwa_mem.py --batchfile DC174_reads.tsv -r Dchry2.2 Dchry2.2.fa --threads 19 --javaXmx 12g --runName default --cleanUpBams --outDir . --tmpDir .
```
