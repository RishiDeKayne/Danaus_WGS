This directory will have all processing scripts for Danaus whole genome paper De-Kayne et al.

All bash scripts can be found in the files [`Danaus_WGS_commands.txt`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/Danaus_WGS_commands.txt) and [`phasing_commands.txt`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/phasing/phasing_commands.txt)

These files explain the following analyses/methods
1. genotype calling
2. PCAs
3. Admixture plots
4. Gemma GWAS
5. Phasing
6. gIMble

<br />  


Additionally there are specific R scripts for the processing, analysing, and plotting output from the follwing steps.
The background file [`DC174_background.csv`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/DC174_background.csv) contains phenotype and indiv. information for all individuals and is used in most R scripts.

**PCA** - [`DC174_PCAs_SVs_only.R`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/DC174_PCAs_SVs_only.R) and [`DC174_makePCAs.R`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/DC174_makePCAs.R)  
**Admixture** - [`DC174_makeAdmixtureplots_sub_regions.R`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/DC174_makeAdmixtureplots_sub_regions.R)  
**GWAS** - Input:[`DC174_make_gemma_input.R`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/DC174_make_gemma_input.R) and Plotting Output: [`DC174_plot_plink_assoc_output.R`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/DC174_plot_plink_assoc_output.R)  
**Population Statistics** - e.g. Input: [`https://github.com/RishiDeKayne/Danaus_WGS/blob/main/make_popsfile.R`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/make_popsfile.R) and Plotting Output: Pi/FST/DXY/Da [`DC174_make_popsfile_plot_pi_50kb.R`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/DC174_make_popsfile_plot_pi_50kb.R)  

<br />  


And various processing scripts including:
[`DC174_get_ploidy.R`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/DC174_get_ploidy.R)  
[`DC174_phase_statistics.R`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/DC174_phase_statistics.R)  
[`DC174_supp_table_background.R`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/DC174_supp_table_background.R)  

<br />  

Finally, the gIMble full model output can be found in [`FINAL_GIMBLE_RESULTS_3morphs.xlsx`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/FINAL_GIMBLE_RESULTS_3morphs.xlsx)

<br />  

LD heatmaps across chr15 were produce using code in [`revision_commands.txt`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/revision_commands.txt)
