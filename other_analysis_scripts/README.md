
### Commands and scripts for various analyses.

Most steps including indel realigtnment, genotype calling, PCA, GWAS, popgen widow statistics and gIMble are detailed in [Danaus_WGS_commands.txt](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/other_analysis_scripts/Danaus_WGS_commands.txt). 

Additionally there are specific R scripts for the processing, analysing, and plotting output from the follwing steps.
The background file [`DC174_background.csv`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/other_analysis_scripts/DC174_background.csv) contains phenotype and indiv. information for all individuals and is used in most R scripts.

**PCA** - [`DC174_PCAs_SVs_only.R`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/other_analysis_scripts/DC174_PCAs_SVs_only.R) and [`DC174_makePCAs.R`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/other_analysis_scripts/DC174_makePCAs.R)  
**Admixture** - [`DC174_makeAdmixtureplots_sub_regions.R`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/other_analysis_scripts/DC174_makeAdmixtureplots_sub_regions.R)  
**GWAS** - Input:[`DC174_make_gemma_input.R`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/other_analysis_scripts/DC174_make_gemma_input.R) and Plotting Output: [`DC174_plot_plink_assoc_output.R`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/other_analysis_scripts/DC174_plot_plink_assoc_output.R)  
**Population Statistics** - e.g. Input: [`make_popsfile.R`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/other_analysis_scripts/make_popsfile.R) and Plotting Output: Pi/FST/DXY/Da [`DC174_make_popsfile_plot_pi_50kb.R`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/other_analysis_scripts/DC174_make_popsfile_plot_pi_50kb.R)  

And various processing scripts including:
[`DC174_get_ploidy.R`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/other_analysis_scripts/DC174_get_ploidy.R)  
[`DC174_phase_statistics.R`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/other_analysis_scripts/DC174_phase_statistics.R)  
[`DC174_supp_table_background.R`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/other_analysis_scripts/DC174_supp_table_background.R)  

Finally, the gIMble full model output can be found in [`FINAL_GIMBLE_RESULTS_3morphs.xlsx`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/other_analysis_scripts/FINAL_GIMBLE_RESULTS_3morphs.xlsx)

LD heatmaps across chr15 were produce using code in [`revision_commands.txt`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/other_analysis_scripts/revision_commands.txt) and plotted in R using [`plot_heat_ma_LD.R`](https://github.com/RishiDeKayne/Danaus_WGS/blob/main/plot_heat_ma_LD.R)
