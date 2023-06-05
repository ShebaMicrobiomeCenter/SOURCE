library(WGCNA)
library(ggplot2)
options(stringsAsFactors = FALSE);

source('/pita/users/tzipi/code/R_figs/WGCNA_funcs.R')

subcohorts_flag = F

min_p = 0.05
min_q = 0.25
# data_type = 'ffq'
# data_type = 'qst'
data_type = 'all'

plot_disease_related_modules_flag = T
# cor_type = 'spearman'
cor_type = 'WGCNA_default'

met_type = 'FFQv10'
# met_type = 'FFQv10Clean'
# met_type = 'SerumMetablolomicsV2LogClean'

in_path = 'rnaSeq_main/res_v2/'
israel_name = sprintf('SOURCE_Israel_TI_with%s_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_%s', met_type, cor_type)
china_name = sprintf('SOURCE_China_TI_ageMatch_with%s_oneSample_pwr_14_netType_signed_hybrid_minMdSz_30_%s', met_type, cor_type)


name = israel_name
israel_ti = load( sprintf('%s/%s/data.RData',in_path, name) )



metadata_prms = c('Gender','Age_years','BMI',
                  'Active_Smoker_manual','Previous_or_Current_Tobacco_Use_manual',
                  'Disease','inflammed','CRP_numeric','Calprotectin_numeric',# 'CRP_log','Calprotectin_log',
                  'health_index_same_sample','health_index_stool')

if (data_type == 'ffq')
{
  moduleTraitCor = moduleTraitCor[, !grepl('^X15',colnames(moduleTraitCor)) ]
  moduleTraitPvalue = moduleTraitPvalue[, !grepl('^X15',colnames(moduleTraitPvalue)) ]
  name = sprintf('%s_%s',name, data_type)
} else if ( data_type == 'qst')
{
  moduleTraitCor = moduleTraitCor[, colnames(moduleTraitCor) %in% metadata_prms | grepl('^X15',colnames(moduleTraitCor)) ]
  moduleTraitPvalue = moduleTraitPvalue[, colnames(moduleTraitPvalue) %in% metadata_prms | grepl('^X15',colnames(moduleTraitPvalue)) ]
  name = sprintf('%s_%s',name, data_type)
}
# out_path = 'rnaSeq_main/res_other_modules/'
# dir.create(out_path)
# out_path = sprintf('%s/china_TI_using_israel_TI_modules_clean_fdr/', out_path)
# dir.create(out_path)



heatmap = WGCNA_ME_met_heatmap(moduleTraitCor, moduleTraitPvalue, name, metadata_prms)
# ggsave(sprintf('%s/%s_module_cor_heatmap.pdf',out_path, name),plot = heatmap, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+dim(moduleTraitPvalue)[2]/10)

## heatmap using predefined module list
# wanted_modules = c('green','yellow','brown','salmon','tan','red','black','purple','pink')
wanted_modules = c('yellow','green','red','pink', 'purple','tan','salmon','black','brown')
wanted_modules = sprintf('ME%s',wanted_modules)
# out_path2 = sprintf('%s/selected_modules/', out_path)
# dir.create(out_path2)
# heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor[wanted_modules,], moduleTraitPvalue[wanted_modules,], name, metadata_prms, wanted_modules_order = wanted_modules )
# ggsave(sprintf('%s/%s_module_cor_heatmap.pdf',out_path2, name),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+dim(moduleTraitPvalue)[2]/10)

# res_p = WGCNA_ME_met_heatmap_full_res(moduleTraitCor[wanted_modules,], moduleTraitPvalue[wanted_modules,], name, 
#                                       metadata_prms, min_p = min_p, wanted_modules_order = wanted_modules, 
#                                       fdr_flag_non_metadata = F, show_p_value_fdr_flag = show_p_value_fdr_flag )
# res_q = WGCNA_ME_met_heatmap_full_res(moduleTraitCor[wanted_modules,], moduleTraitPvalue[wanted_modules,], name, 
#                                       metadata_prms, min_p = min_q, wanted_modules_order = wanted_modules, 
#                                       fdr_flag_non_metadata = T, show_p_value_fdr_flag = show_p_value_fdr_flag )
res = WGCNA_ME_met_heatmap_pq(moduleTraitCor, moduleTraitPvalue, name, metadata_prms = metadata_prms, 
                                   min_p = min_p, wanted_modules_order = wanted_modules, 
                                   cluster_cols_order_flag = T, show_p_value_fdr_flag = T, min_q = min_q)
heatmap_s_israel = res[[1]]
y_lables = ggplot_build(heatmap_s_israel)$layout$panel_params[[1]]$y$get_labels()

out_path2 = sprintf('%s/%s/selected_modules',in_path, name )
dir.create(out_path2)
ggsave(sprintf('%s/%s_module_cor_heatmap_p%s_q%s.pdf',out_path2, name, min_p, min_q),plot = heatmap_s_israel, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+length(y_lables)/10)

mdc = res[[2]]
write.table(x = mdc, file = sprintf('%s/%s_module_cor_df.pdf',out_path2, name), quote = F, sep = '\t', row.names = F)


## save israel data
israel_df = df
israel_name = name


## china side
# name = 'SOURCE_China_TI_withFFQKcalNorm_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_spearman'
# name = 'SOURCE_China_TI_withFFQabsolute_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_spearman'
name = china_name
china_ti = load( sprintf('%s/%s/data.RData',in_path, name) )


## sync genes between the datasets
good_genes = intersect(names(ftr_df), israel_df$Gene)
israel_df = israel_df[israel_df$Gene %in% good_genes, ]
ftr_df = ftr_df[, israel_df$Gene]

## ME calculation
# Calculate MEs with color labels
# name = sprintf('SOURCE_China_TI_by_%s', israel_name)
name = sprintf('SOURCE_China_ageMatch_TI_by_%s', israel_name)
# WGCNA_plot_genes_correlaitons_by_module(ftr_df, israel_df$module_color, out_path)
res = WGCNA_calculate_module_eigengenes_and_corrs(ftr_df,israel_df$module_color, metadata_df, cor_type, name)
MEs = res[[1]]
moduleTraitCor = res[[2]]
moduleTraitPvalue = res[[3]]

if (data_type == 'ffq')
{
  moduleTraitCor = moduleTraitCor[, !grepl('^X15',colnames(moduleTraitCor)) ]
  moduleTraitPvalue = moduleTraitPvalue[, !grepl('^X15',colnames(moduleTraitPvalue)) ]
} else if ( data_type == 'qst')
{
  moduleTraitCor = moduleTraitCor[, colnames(moduleTraitCor) %in% metadata_prms | grepl('^X15',colnames(moduleTraitCor)) ]
  moduleTraitPvalue = moduleTraitPvalue[, colnames(moduleTraitPvalue) %in% metadata_prms | grepl('^X15',colnames(moduleTraitPvalue)) ]
}

heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor[wanted_modules,], moduleTraitPvalue[wanted_modules,], name, metadata_prms, min_p = min_p, wanted_modules_order = wanted_modules, fdr_flag_non_metadata = F,  show_p_value_fdr_flag = T)
y_lables = ggplot_build(heatmap_s)$layout$panel_params[[1]]$y$get_labels()
# ggsave(sprintf('%s/%s_module_cor_heatmap_p%s.pdf',out_path2, name, min_p),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+length(y_lables)/10)

res = WGCNA_ME_met_heatmap_pq(moduleTraitCor[wanted_modules,], moduleTraitPvalue[wanted_modules,], name, metadata_prms = metadata_prms,
                              min_p = min_p, wanted_modules_order = wanted_modules,
                              cluster_cols_order_flag = T, show_p_value_fdr_flag = T, min_q = min_q)
heatmap_s_china = res[[1]]


heatmap_s_israel
heatmap_s_china





