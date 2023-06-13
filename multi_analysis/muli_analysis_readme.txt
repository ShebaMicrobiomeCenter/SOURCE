
This is a guide for the multi_analysis section of this code.
Here multiple omics levels data is integrated, using WGCNA and halla as the base analysis. Details abour intiall analysis within each omic level are under the appropriate subdirectory. 
The relevant neccesary functions are under ../funcs. hard paths were not all rewirten in the code. 

## WGCNA (Weighted Correlation Network Analysis)

RNAseq data was used as the base for WCNA analysis, and later correalted to diet and metabolimics data.
First run WGCNA/WGCNA_base_run.R scripts for both Israel and China data. 
Than use WGCNA/WGCNA_using_other_modules_fdr.R to apply israel modules to China data. 
The rest of the dowstream analysis is uder WGCNA/ with standalone meaningfully named scripts.  
functionla enrichment of the genes associated with each module was performed using toppgene (online), the code for plotting this enrichment is under WGCNA/modules_enrichment/


## HAllA (High-sensitivity pattern discovery in large, paired multi-omic datasets)

Omics pairs correlations were analysed usign HAllA. 
Initial set-up of paired data was performed using the scripts:
set_halla_input_prepairs.R
fit_halla_input_pairs_v2.R

and actual HAllA run with the script:
halla_script.sh

The metabolites used for some of the HAllA aalysis are the metabolites associated with Israel ileal transcriptomics WGCNA modules, calculated in the WGCNA section of the code. 

The rest of the dowstream analysis is uder halla/ with standalone meaningfully named scripts. 


