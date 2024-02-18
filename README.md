# SOURCE

This is the code associated with the paper : 

Diet-omics in the Study of Urban and Rural Crohn Disease Evolution (SOURCE)


Crohn disease (CD) burden has increased with globalization/urbanization, and with the rapid pace of the increase is attributed to environmental changes rather than genetic drift. The Study Of Urban and Rural CD Evolution (SOURCE, n=380) considered diet-omics domains simultaneously to detect complex interactions and identify potential beneficial and pathogenic factors linked with rural-urban transition and CD. We characterized exposures, diet, ileal transcriptomics, metabolomics, and microbiome in newly diagnosed Crohn Disease (CD) patients and controls in rural and urban China and Israel. Time spent by rural residents in urban environments was linked with changes in gut microbial composition and metabolomics, which mirrored those seen in CD. Ileal transcriptomics highlighted host metabolic and immune gene expression modules, that were linked to potential protective dietary exposures (coffee, manganese, vitamin D), fecal metabolites, and the microbiome. Bacteria-associated metabolites were primarily linked with host immune modules, whereas diet-linked metabolites were associated with host epithelial metabolic functions. 


The functions used through the code should all be available under the funcs/ directory (with some additional functions that were not eventually used here).  
Example data files are under the data/ subdirectory within the appropriate directory.


Additional run guides are available as readme.txt files under the appropriate subdirectories, as necessary. For example, the 16S run guide is under 16S/16S_readme.txt

For the code of the multi-omic analysis performed by the Borenstein lab, please visit their [GitHub Repository](https://github.com/borenstein-lab/SOURCE_multiomics).


### Raw data is publicly available at the following:
RNASeq Israel and China datasets were deposited in GEO: GSE199906 and GSE233900, respectively.

The 16S amplicon sequencing dataset was deposited in the National Center for Biotechnology Information as BioProject PRJNA978342.


RNASeq Israel and China datasets generated in this study have been deposited in the GEO database under accession code: [GSE199906](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199906) and [GSE233900](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233900), respectively.

The 16S amplicon sequencing dataset generated in this study has been deposited in BioProject under accession code: [PRJNA978342](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA978342).

Microbial shotgun sequencing generated in this study have been deposited BioProject under accession code: [PRJNA1056458](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1056458).



## R session info:
```
R version 4.1.2 (2021-11-01)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.3 LTS

Matrix products: default
BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggraph_2.1.0                ggalluvial_0.12.5           networkD3_0.4              
 [4] RColorBrewer_1.1-3          DESeq2_1.34.0               SummarizedExperiment_1.24.0
 [7] Biobase_2.54.0              MatrixGenerics_1.6.0        GenomicRanges_1.46.1       
[10] GenomeInfoDb_1.30.1         IRanges_2.28.0              S4Vectors_0.32.4           
[13] BiocGenerics_0.40.0         tximport_1.22.0             ggmosaic_0.3.3             
[16] Hmisc_4.7-2                 Formula_1.2-4               survival_3.2-13            
[19] ggside_0.2.2                Maaslin2_1.8.0              vegan_2.6-4                
[22] lattice_0.20-45             permute_0.9-7               matrixStats_0.63.0         
[25] stringr_1.5.0               reshape_0.8.9               patchwork_1.1.2            
[28] ggplot2_3.4.1               WGCNA_1.72-1                fastcluster_1.2.3          
[31] dynamicTreeCut_1.63-1      

loaded via a namespace (and not attached):
  [1] colorspace_2.1-0       ggsignif_0.6.4         deldir_1.0-6          
  [4] htmlTable_2.4.1        XVector_0.34.0         base64enc_0.1-3       
  [7] rstudioapi_0.14        farver_2.1.1           graphlayouts_0.8.4    
 [10] getopt_1.20.3          ggrepel_0.9.2          bit64_4.0.5           
 [13] AnnotationDbi_1.56.2   fansi_1.0.4            mvtnorm_1.1-3         
 [16] codetools_0.2-18       splines_4.1.2          doParallel_1.0.17     
 [19] cachem_1.0.6           impute_1.68.0          robustbase_0.95-0     
 [22] geneplotter_1.72.0     knitr_1.41             polyclip_1.10-4       
 [25] jsonlite_1.8.4         annotate_1.72.0        cluster_2.1.2         
 [28] GO.db_3.14.0           png_0.1-8              ggforce_0.4.1         
 [31] compiler_4.1.2         httr_1.4.4             backports_1.4.1       
 [34] Matrix_1.5-3           fastmap_1.1.0          lazyeval_0.2.2        
 [37] cli_3.6.0              tweenr_2.0.2           htmltools_0.5.4       
 [40] tools_4.1.2            igraph_1.4.1           gtable_0.3.1          
 [43] glue_1.6.2             GenomeInfoDbData_1.2.7 reshape2_1.4.4        
 [46] dplyr_1.1.0            Rcpp_1.0.10            biglm_0.9-2.1         
 [49] vctrs_0.5.2            Biostrings_2.62.0      preprocessCore_1.56.0 
 [52] nlme_3.1-153           iterators_1.0.14       optparse_1.7.3        
 [55] xfun_0.36              lifecycle_1.0.3        XML_3.99-0.13         
 [58] DEoptimR_1.0-11        zlibbioc_1.40.0        MASS_7.3-54           
 [61] scales_1.2.1           tidygraph_1.2.3        parallel_4.1.2        
 [64] memoise_2.0.1          gridExtra_2.3          rpart_4.1-15          
 [67] latticeExtra_0.6-30    stringi_1.7.12         RSQLite_2.2.20        
 [70] genefilter_1.76.0      pcaPP_2.0-3            foreach_1.5.2         
 [73] checkmate_2.1.0        BiocParallel_1.28.3    rlang_1.0.6           
 [76] pkgconfig_2.0.3        bitops_1.0-7           lpsymphony_1.22.0     
 [79] purrr_1.0.1            htmlwidgets_1.6.1      bit_4.0.5             
 [82] tidyselect_1.2.0       plyr_1.8.8             magrittr_2.0.3        
 [85] R6_2.5.1               generics_0.1.3         DelayedArray_0.20.0   
 [88] DBI_1.1.3              pillar_1.8.1           foreign_0.8-81        
 [91] withr_2.5.0            mgcv_1.8-38            KEGGREST_1.34.0       
 [94] RCurl_1.98-1.9         nnet_7.3-16            tibble_3.1.8          
 [97] crayon_1.5.2           interp_1.1-3           utf8_1.2.3            
[100] plotly_4.10.1          viridis_0.6.2          jpeg_0.1-10           
[103] locfit_1.5-9.7         grid_4.1.2             data.table_1.14.6     
[106] blob_1.2.3             digest_0.6.31          xtable_1.8-4          
[109] tidyr_1.3.0            munsell_0.5.0          viridisLite_0.4.1 
```
