library(WGCNA)
library(ggplot2)
options(stringsAsFactors = FALSE);

source('/pita/users/tzipi/code/R_figs/WGCNA_funcs.R')

in_path = 'rnaSeq_main/res_v2/'
# name = 'SOURCE_Israel_TI_withStoolMetablolomics_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default'
# name = 'SOURCE_Israel_R_withFFQv10_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default'
# name = 'SOURCE_Israel_TI_withStoolMetablolomicsLogV2_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default'
name = 'SOURCE_Israel_TI_withStoolMetablolomicsV2LogClean_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default'
name = 'SOURCE_Israel_TI_withStoolMetablolomicsV2Log_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default'
name = 'SOURCE_Israel_TI_withStool16S_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default'
rdata_names = load( sprintf('%s/%s/data.RData',in_path, name) )


 
check_ftr = 'Cytosine_112.1'
check_module = 'yellow'

check_ftr = 'f__Bifidobacteriaceae.g__Bifidobacterium_ASV11475'
check_module = 'salmon'



check_module = sprintf('ME%s', check_module)

if ( ! all(row.names(MEs) == row.names(metadata_df)) )
  print('Data sync problem!')

# cor_res = cor.test(x = metadata_df[[check_ftr]], y=MEs[[check_module]], method = 'pearson')
cor_res = cor.test(x = metadata_df[[check_ftr]], y=MEs[[check_module]], method = 'spearman')
g = ggplot(metadata_df, aes(x=metadata_df[[check_ftr]], y= MEs[[check_module]], colour = as.factor(Disease))) +
  geom_point() + xlab(check_ftr) + ylab(check_module) + theme_bw()+  
  ggtitle(sprintf('%s\nr = %.3f, p = %s', cor_res$method, cor_res$estimate, signif(cor_res$p.value, 1)))


# cor_res = cor.test(x = metadata_df$`4-Hydroxy-L-Phenylglycine_168.1`, y=MEs[[check_module]], method = 'pearson')
# g = ggplot(metadata_df, aes(x=metadata_df$`4-Hydroxy-L-Phenylglycine_168.1`, y= MEs[[check_module]], colour = as.factor(Disease))) +
#   geom_point() + xlab(check_ftr) + ylab(check_module) + theme_bw()+  
#   ggtitle(sprintf('%s\nr = %.3f, p = %s', cor_res$method, cor_res$estimate, signif(cor_res$p.value, 1)))



