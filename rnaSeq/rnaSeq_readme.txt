
This is a guide for the RNAseq section of this code.
RNAseq intergration with other omics is found under multi_analysis. 
The relevant neccesary functions are all under ../funcs. hard paths were not all rewirten in the code. 

For the initial kallisto run starting from raw fastq files use the script: 
run_kallisto.py.

transcripts are than summarized to gene level using: 
trans2gene_tximport.R. 

PCA plots are gnerated using:
pca_v2.R
