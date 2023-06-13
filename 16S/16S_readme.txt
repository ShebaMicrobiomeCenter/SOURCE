
This is a guide for the 16S section of this code.
16S intergration with other omics is found under multi_analysis. Here you can find the code only for 16S with FFQ+questioneaire data analysis. 
The relevant neccesary functions are all under ../funcs. hard paths were not all rewirten in the code. 

## QIIME2 initial analysis

QIIME2 analysis was all performed using version 2021.4.

For the initial QIIME2 run starting from raw fastq files use the scripts:
set_data/qiime2_basic_run.py

biom is than filtered to remove contaminents using calour/dbbact (code currently missing, will be updated), and then taxonomy and diversity analysis id performed using the scripts:

set_data/run_amnon_filtered_biom_for_permanova.py
set_data/beta_dist_qiime2_script.sh


## R downstream analysis

first, the following scripts are run to set data tables in a format easy to use in dowstream analysis:
R/update_vars_qiime2.R
R/update_vars_qiime2_SVs.R 

The rest of the dowstream analysis is also uder R/ with standalone meaningfully named scripts. The exceptiosn are:
plot_maaslin2.R and fig_rural_maas_res.R, that requires maaslin2.R to be run first. 
plot_16S_permanova.R, that requires permonova_16S_v2.R be run first. 
 
