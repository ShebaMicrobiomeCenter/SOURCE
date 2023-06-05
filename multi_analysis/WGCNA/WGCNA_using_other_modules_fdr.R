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

in_path = 'rnaSeq_main/res_v2/'
out_path = 'rnaSeq_main/res_other_modules/'
dir.create(out_path)

# out_path = sprintf('%s/TI_china_israel_modules_cleanFFQ/', out_path)
# israel_name = sprintf('SOURCE_Israel_TI_withFFQv10Clean_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_%s',cor_type)
# china_name = sprintf('SOURCE_China_TI_ageMatch_withFFQv10Clean_oneSample_pwr_14_netType_signed_hybrid_minMdSz_30_%s',cor_type)

out_path = sprintf('%s/TI_china_israel_modules_FFQ/', out_path)
israel_name = sprintf('SOURCE_Israel_TI_withFFQv10_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_%s',cor_type)
china_name = sprintf('SOURCE_China_TI_ageMatch_withFFQv10_oneSample_pwr_14_netType_signed_hybrid_minMdSz_30_%s',cor_type)

dir.create(out_path)

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



heatmap = WGCNA_ME_met_heatmap(moduleTraitCor, moduleTraitPvalue, name, metadata_prms)
ggsave(sprintf('%s/%s_module_cor_heatmap.pdf',out_path, name),plot = heatmap, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+dim(moduleTraitPvalue)[2]/10)

## heatmap using predefined module list
# wanted_modules = c('green','yellow','brown','salmon','tan','red','black','purple','pink')
wanted_modules = c('yellow','green','red','pink', 'purple','tan','salmon','black','brown')
wanted_modules = sprintf('ME%s',wanted_modules)
out_path2 = sprintf('%s/selected_modules/', out_path)
dir.create(out_path2)
heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor[wanted_modules,], moduleTraitPvalue[wanted_modules,], name, metadata_prms, wanted_modules_order = wanted_modules )
ggsave(sprintf('%s/%s_module_cor_heatmap.pdf',out_path2, name),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+dim(moduleTraitPvalue)[2]/10)

# res = WGCNA_ME_met_heatmap_full_res(moduleTraitCor[wanted_modules,], moduleTraitPvalue[wanted_modules,], name, metadata_prms, min_p = min_q, wanted_modules_order = wanted_modules, fdr_flag_non_metadata = T, show_p_value_fdr_flag = T )
# heatmap_s = res[[1]]
# y_lables = ggplot_build(heatmap_s)$layout$panel_params[[1]]$y$get_labels()
# ggsave(sprintf('%s/%s_module_cor_heatmap_q%s.pdf',out_path2, name, min_q),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+length(y_lables)/10)

# res = WGCNA_ME_met_heatmap_full_res(moduleTraitCor[wanted_modules,], moduleTraitPvalue[wanted_modules,], name, metadata_prms, min_p = min_p, wanted_modules_order = wanted_modules, fdr_flag_non_metadata = F, show_p_value_fdr_flag = T )
# heatmap_s = res[[1]]
# y_lables = ggplot_build(heatmap_s)$layout$panel_params[[1]]$y$get_labels()
# ggsave(sprintf('%s/%s_module_cor_heatmap_p%s.pdf',out_path2, name, min_p),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+length(y_lables)/10)

res_israel = WGCNA_ME_met_heatmap_pq(moduleTraitCor[wanted_modules,], moduleTraitPvalue[wanted_modules,], name, metadata_prms = metadata_prms, 
                              min_p = min_p, wanted_modules_order = wanted_modules, 
                              cluster_cols_order_flag = T, show_p_value_fdr_flag = T, min_q = min_q)
heatmap_s = res_israel[[1]]
y_lables = ggplot_build(heatmap_s)$layout$panel_params[[1]]$y$get_labels()
ggsave(sprintf('%s/%s_module_cor_heatmap_p%s_q%s.pdf',out_path2, name, min_p, min_q),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+length(y_lables)/10)

if (subcohorts_flag)
{
  ## for looking at only a subcohort ( only diease/healthy)
  subcohort_name = 'Disease'
  wanted_pos = metadata_df$Disease == 1
  
  res = WGCNA_calculate_module_eigengenes_and_corrs(ftr_df[wanted_pos,],moduleColors, clean_metadata_eampty_noVar( metadata_df[wanted_pos,] ), cor_type, name)
  moduleTraitCor_d = res[[2]]
  moduleTraitPvalue_d = res[[3]]
  wanted_labels = y_lables[y_lables %in% colnames(moduleTraitCor_d)]
  
  heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor_d[wanted_modules,], moduleTraitPvalue_d[wanted_modules,], name, metadata_prms, wanted_modules_order = wanted_modules)
  ggsave(sprintf('%s/%s_module_cor_heatmap_%s.pdf',out_path2, name,subcohort_name),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue_d)[1]/2,height = 5+dim(moduleTraitPvalue_d)[2]/10)
  
  heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor_d[wanted_modules,wanted_labels], moduleTraitPvalue_d[wanted_modules,wanted_labels], name, metadata_prms, wanted_modules_order = wanted_modules, cluster_cols_order_flag = F )
  ggsave(sprintf('%s/%s_module_cor_heatmap_%s_allq%s.pdf',out_path2, name,subcohort_name, min_q),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue_d)[1]/2,height = 5+length(wanted_labels)/10)
  
  
  subcohort_name = 'Healthy'
  wanted_pos = metadata_df$Disease == 0
  
  res = WGCNA_calculate_module_eigengenes_and_corrs(ftr_df[wanted_pos,],moduleColors, clean_metadata_eampty_noVar( metadata_df[wanted_pos,] ), cor_type, name)
  moduleTraitCor_d = res[[2]]
  moduleTraitPvalue_d = res[[3]]
  wanted_labels = y_lables[y_lables %in% colnames(moduleTraitCor_d)]
  
  heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor_d[wanted_modules,], moduleTraitPvalue_d[wanted_modules,], name, metadata_prms, wanted_modules_order = wanted_modules)
  ggsave(sprintf('%s/%s_module_cor_heatmap_%s.pdf',out_path2, name,subcohort_name),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue_d)[1]/2,height = 5+dim(moduleTraitPvalue_d)[2]/10)
  
  heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor_d[wanted_modules,wanted_labels], moduleTraitPvalue_d[wanted_modules,wanted_labels], name, metadata_prms, wanted_modules_order = wanted_modules, cluster_cols_order_flag = F )
  ggsave(sprintf('%s/%s_module_cor_heatmap_%s_allq%s.pdf',out_path2, name,subcohort_name, min_q),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue_d)[1]/2,height = 5+length(wanted_labels)/10)
  
}

## save israel data
israel_df = df
# israel_name = name


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
  name = sprintf('%s_%s',name, data_type)
} else if ( data_type == 'qst')
{
  moduleTraitCor = moduleTraitCor[, colnames(moduleTraitCor) %in% metadata_prms | grepl('^X15',colnames(moduleTraitCor)) ]
  moduleTraitPvalue = moduleTraitPvalue[, colnames(moduleTraitPvalue) %in% metadata_prms | grepl('^X15',colnames(moduleTraitPvalue)) ]
  name = sprintf('%s_%s',name, data_type)
}

# my version
heatmap = WGCNA_ME_met_heatmap(moduleTraitCor, moduleTraitPvalue, name, metadata_prms)
ggsave(sprintf('%s/%s_module_cor_heatmap.pdf',out_path, name),plot = heatmap, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+dim(moduleTraitPvalue)[2]/10)
if (plot_disease_related_modules_flag)
{
  heatmap2 = WGCNA_ME_met_heatmap_prm_related(moduleTraitCor, moduleTraitPvalue, name, prm = 'Disease', cutoff = min_p, metadata_prms)
  ggsave(sprintf('%s/%s_module_cor_heatmap_%sF%s.pdf',out_path, name,'Disease',min_p), plot = heatmap2, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+dim(moduleTraitPvalue)[2]/10)
}


# heatmap using predefined module list
heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor[wanted_modules,], moduleTraitPvalue[wanted_modules,], name, metadata_prms, wanted_modules_order = wanted_modules)
ggsave(sprintf('%s/%s_module_cor_heatmap.pdf',out_path2, name),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+dim(moduleTraitPvalue)[2]/10)

# heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor[wanted_modules,], moduleTraitPvalue[wanted_modules,], name, metadata_prms, min_p = min_q, wanted_modules_order = wanted_modules, fdr_flag_non_metadata = T,  show_p_value_fdr_flag = T)
# y_lables = ggplot_build(heatmap_s)$layout$panel_params[[1]]$y$get_labels()
# ggsave(sprintf('%s/%s_module_cor_heatmap_q%s.pdf',out_path2, name, min_q),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+length(y_lables)/10)

res_china = WGCNA_ME_met_heatmap_pq(moduleTraitCor[wanted_modules,], moduleTraitPvalue[wanted_modules,], name, metadata_prms = metadata_prms, 
                              min_p = min_p, wanted_modules_order = wanted_modules, 
                              cluster_cols_order_flag = T, show_p_value_fdr_flag = T, min_q = min_q)
heatmap_s = res_china[[1]]
y_lables = ggplot_build(heatmap_s)$layout$panel_params[[1]]$y$get_labels()
ggsave(sprintf('%s/%s_module_cor_heatmap_p%s_q%s.pdf',out_path2, name, min_p, min_q),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+length(y_lables)/10)

# heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor[wanted_modules,], moduleTraitPvalue[wanted_modules,], name, metadata_prms, min_p = min_p, wanted_modules_order = wanted_modules, fdr_flag_non_metadata = F,  show_p_value_fdr_flag = T)
# y_lables = ggplot_build(heatmap_s)$layout$panel_params[[1]]$y$get_labels()
# ggsave(sprintf('%s/%s_module_cor_heatmap_p%s.pdf',out_path2, name, min_p),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+length(y_lables)/10)

if (subcohorts_flag)
{
  # for looking at only a subcohort ( only diease/healthy)
  subcohort_name = 'Disease'
  wanted_pos = metadata_df$Disease == 1
  
  res = WGCNA_calculate_module_eigengenes_and_corrs(ftr_df[wanted_pos,],israel_df$module_color, clean_metadata_eampty_noVar( metadata_df[wanted_pos,] ), cor_type, name)
  moduleTraitCor_d = res[[2]]
  moduleTraitPvalue_d = res[[3]]
  wanted_labels = y_lables[y_lables %in% colnames(moduleTraitCor_d)]
  
  heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor_d[wanted_modules,], moduleTraitPvalue_d[wanted_modules,], name, metadata_prms, wanted_modules_order = wanted_modules)
  ggsave(sprintf('%s/%s_module_cor_heatmap_%s.pdf',out_path2, name,subcohort_name),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue_d)[1]/2,height = 5+dim(moduleTraitPvalue_d)[2]/10)
  
  heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor_d[wanted_modules,wanted_labels], moduleTraitPvalue_d[wanted_modules,wanted_labels], name, metadata_prms, wanted_modules_order = wanted_modules, cluster_cols_order_flag = F )
  ggsave(sprintf('%s/%s_module_cor_heatmap_%s_allq%s.pdf',out_path2, name,subcohort_name, min_q),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue_d)[1]/2,height = 5+length(wanted_labels)/10)
  
  
  subcohort_name = 'Healthy'
  wanted_pos = metadata_df$Disease == 0
  
  res = WGCNA_calculate_module_eigengenes_and_corrs(ftr_df[wanted_pos,],israel_df$module_color, clean_metadata_eampty_noVar( metadata_df[wanted_pos,] ), cor_type, name)
  moduleTraitCor_d = res[[2]]
  moduleTraitPvalue_d = res[[3]]
  wanted_labels = y_lables[y_lables %in% colnames(moduleTraitCor_d)]
  
  heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor_d[wanted_modules,], moduleTraitPvalue_d[wanted_modules,], name, metadata_prms, wanted_modules_order = wanted_modules)
  # ggsave(sprintf('%s/%s_module_cor_heatmap_%s.pdf',out_path2, name,subcohort_name),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue_d)[1]/2,height = 5+dim(moduleTraitPvalue_d)[2]/10)
  
  heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor_d[wanted_modules,wanted_labels], moduleTraitPvalue_d[wanted_modules,wanted_labels], name, metadata_prms, wanted_modules_order = wanted_modules, cluster_cols_order_flag = F )
  # ggsave(sprintf('%s/%s_module_cor_heatmap_%s_allq%s.pdf',out_path2, name,subcohort_name, min_q),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue_d)[1]/2,height = 5+length(wanted_labels)/10)
  
  
  # g = ggplot(metadata_df, aes(x=as.factor(metadata_df$Active_Smoker_manual), y=MEs$MEtan)) +
  #   geom_boxplot(outlier.alpha = 0) + geom_jitter(aes(colour = as.factor(metadata_df$Disease))) +
  #   scale_color_manual(values = c('blue','red'), name = 'Disease') + xlab('Active_Smoker_manual') + # ylab('tan') +
  #   theme_bw()
  # 
  # ggsave(sprintf('%s/%s_smoking_active_tan_box.pdf',out_path2, name),plot = g, device = 'pdf', width = 4,height = 3.5)
}
