
This is a guide for the metabolomics section of this code.
metabolomics intergration with other omics is found under multi_analysis. 
The relevant neccesary functions are all under ../funcs. hard paths were not all rewirten in the code. 

## initial analysis

basic analysis was preformed using MZmine2.53 for the Israeli data, and MassLynx software v4.1 for the China data, no code available.
data was than normalized using norm_data.R 


## downstream analysis

PCoA canberra analysis was initally performed using qiime-2021.4, with the script set_data_qiime2.sh
final ploting and statistics was than performed in R with plot_PCoA.R, plot_PCoA_israel.R

maaslin2 basic run was performed with masslin2.R, further analysis and ploting was performed with the scripts fig_rural_maas_res.R, sactter_rural_dx_maas_coef.R

