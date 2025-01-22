
### Computation of the Neighbor-Net network. 

Distance matrix was produced using the script `distMat.py` available here: https://github.com/simonhmartin/genomics_general

```bash
tabix -h DC174OG.contig15.1.realigned.DP7GQ30het75min87minVar2.geno.gz contig15.1:4035981-5322256 contig15.1:5322257-6220875 contig15.1:6251636-7829094 contig15.1:14803486-15299732 > temp.geno && python ~/Research/genomics_general/distMat.py -g temp.geno -o DC174OG.contig15.1.realigned.DP7GQ30het75min87minVar2.regions1-4.dist -f phased --windType cat && rm temp.geno
```

The distance matrix was then loaded into SplitsTree manually to make the svg.
