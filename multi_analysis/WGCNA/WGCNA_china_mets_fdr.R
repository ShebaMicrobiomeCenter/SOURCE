library(WGCNA)
library(ggplot2)
library(patchwork)
options(stringsAsFactors = FALSE);

source('/pita/users/tzipi/code/R_figs/WGCNA_funcs.R')

subcohorts_flag = T

min_p = 0.05
min_q = 0.25
# data_type = 'ffq'
# data_type = 'qst'
data_type = 'all'

col2 = colorRampPalette(rev( c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                               "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061") ))(200)

plot_disease_related_modules_flag = T
# cor_type = 'spearman'
cor_type = 'WGCNA_default'

in_path = 'rnaSeq_main/res_v2/'
out_path = 'rnaSeq_main/res_other_modules/'
dir.create(out_path)

# out_path = sprintf('%s/TI_china_israel_modules_cleanFFQ/', out_path)
# israel_name = sprintf('SOURCE_Israel_TI_withFFQv10Clean_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_%s',cor_type)
# china_name = sprintf('SOURCE_China_TI_ageMatch_withFFQv10Clean_oneSample_pwr_14_netType_signed_hybrid_minMdSz_30_%s',cor_type)

out_path = sprintf('%s/TI_china_israel_modules_metsLogClean/', out_path)
israel_name = sprintf('SOURCE_Israel_TI_withFFQv10_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_%s',cor_type)
china_name = sprintf('SOURCE_China_TI_ageMatch_withFFQv10_oneSample_pwr_14_netType_signed_hybrid_minMdSz_30_%s',cor_type)
wanted_modules = c('yellow','green','red','pink', 'purple','tan','salmon','black','brown')
tt='TI'

# out_path = sprintf('%s/R_china_israel_modules_metsLogClean/', out_path)
# israel_name = sprintf('SOURCE_Israel_R_withFFQv10_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_%s',cor_type)
# china_name = sprintf('SOURCE_China_R_ageMatch_withFFQv10_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_%s',cor_type)
# wanted_modules = c('brown','tan','purple','greenyellow','red','turquoise','black','green')
# tt = 'R'

dir.create(out_path)

name = israel_name

israel_ti = load( sprintf('%s/%s/data.RData',in_path, name) )

metadata_prms = c('Gender','Age_years','BMI',
                  'Active_Smoker_manual','Previous_or_Current_Tobacco_Use_manual',
                  'Disease','inflammed','CRP_numeric','Calprotectin_numeric',# 'CRP_log','Calprotectin_log',
                  'health_index_same_sample','health_index_stool')

# heatmap = WGCNA_ME_met_heatmap(moduleTraitCor, moduleTraitPvalue, name, metadata_prms)
# ggsave(sprintf('%s/%s_module_cor_heatmap.pdf',out_path, name),plot = heatmap, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+dim(moduleTraitPvalue)[2]/10)

## heatmap using predefined module list
# wanted_modules = c('green','yellow','brown','salmon','tan','red','black','purple','pink')
wanted_modules = sprintf('ME%s',wanted_modules)
out_path2 = sprintf('%s/selected_modules/', out_path)
dir.create(out_path2)
# heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor[wanted_modules,], moduleTraitPvalue[wanted_modules,], name, metadata_prms, wanted_modules_order = wanted_modules )
# ggsave(sprintf('%s/%s_module_cor_heatmap.pdf',out_path2, name),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+dim(moduleTraitPvalue)[2]/10)

# res = WGCNA_ME_met_heatmap_full_res(moduleTraitCor[wanted_modules,], moduleTraitPvalue[wanted_modules,], name, metadata_prms, min_p = min_q, wanted_modules_order = wanted_modules, fdr_flag_non_metadata = T, show_p_value_fdr_flag = T )
# heatmap_s = res[[1]]
# y_lables = ggplot_build(heatmap_s)$layout$panel_params[[1]]$y$get_labels()
# ggsave(sprintf('%s/%s_module_cor_heatmap_q%s.pdf',out_path2, name, min_q),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+length(y_lables)/10)

# res = WGCNA_ME_met_heatmap_full_res(moduleTraitCor[wanted_modules,], moduleTraitPvalue[wanted_modules,], name, metadata_prms, min_p = min_p, wanted_modules_order = wanted_modules, fdr_flag_non_metadata = F, show_p_value_fdr_flag = T )
# heatmap_s = res[[1]]
# y_lables = ggplot_build(heatmap_s)$layout$panel_params[[1]]$y$get_labels()
# ggsave(sprintf('%s/%s_module_cor_heatmap_p%s.pdf',out_path2, name, min_p),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+length(y_lables)/10)

# res_israel = WGCNA_ME_met_heatmap_pq(moduleTraitCor[wanted_modules,], moduleTraitPvalue[wanted_modules,], name, metadata_prms = metadata_prms, 
#                                      min_p = min_p, wanted_modules_order = wanted_modules, 
#                                      cluster_cols_order_flag = T, show_p_value_fdr_flag = T, min_q = min_q)
# heatmap_s = res_israel[[1]]
# y_lables = ggplot_build(heatmap_s)$layout$panel_params[[1]]$y$get_labels()
# ggsave(sprintf('%s/%s_module_cor_heatmap_p%s_q%s.pdf',out_path2, name, min_p, min_q),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+length(y_lables)/10)

# if (subcohorts_flag)
# {
#   ## for looking at only a subcohort ( only diease/healthy)
#   subcohort_name = 'Disease'
#   wanted_pos = metadata_df$Disease == 1
#   
#   res = WGCNA_calculate_module_eigengenes_and_corrs(ftr_df[wanted_pos,],moduleColors, clean_metadata_eampty_noVar( metadata_df[wanted_pos,] ), cor_type, name)
#   moduleTraitCor_d = res[[2]]
#   moduleTraitPvalue_d = res[[3]]
#   wanted_labels = y_lables[y_lables %in% colnames(moduleTraitCor_d)]
#   
#   heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor_d[wanted_modules,], moduleTraitPvalue_d[wanted_modules,], name, metadata_prms, wanted_modules_order = wanted_modules)
#   ggsave(sprintf('%s/%s_module_cor_heatmap_%s.pdf',out_path2, name,subcohort_name),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue_d)[1]/2,height = 5+dim(moduleTraitPvalue_d)[2]/10)
#   
#   heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor_d[wanted_modules,wanted_labels], moduleTraitPvalue_d[wanted_modules,wanted_labels], name, metadata_prms, wanted_modules_order = wanted_modules, cluster_cols_order_flag = F )
#   ggsave(sprintf('%s/%s_module_cor_heatmap_%s_allq%s.pdf',out_path2, name,subcohort_name, min_q),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue_d)[1]/2,height = 5+length(wanted_labels)/10)
#   
#   
#   subcohort_name = 'Healthy'
#   wanted_pos = metadata_df$Disease == 0
#   
#   res = WGCNA_calculate_module_eigengenes_and_corrs(ftr_df[wanted_pos,],moduleColors, clean_metadata_eampty_noVar( metadata_df[wanted_pos,] ), cor_type, name)
#   moduleTraitCor_d = res[[2]]
#   moduleTraitPvalue_d = res[[3]]
#   wanted_labels = y_lables[y_lables %in% colnames(moduleTraitCor_d)]
#   
#   heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor_d[wanted_modules,], moduleTraitPvalue_d[wanted_modules,], name, metadata_prms, wanted_modules_order = wanted_modules)
#   ggsave(sprintf('%s/%s_module_cor_heatmap_%s.pdf',out_path2, name,subcohort_name),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue_d)[1]/2,height = 5+dim(moduleTraitPvalue_d)[2]/10)
#   
#   heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor_d[wanted_modules,wanted_labels], moduleTraitPvalue_d[wanted_modules,wanted_labels], name, metadata_prms, wanted_modules_order = wanted_modules, cluster_cols_order_flag = F )
#   ggsave(sprintf('%s/%s_module_cor_heatmap_%s_allq%s.pdf',out_path2, name,subcohort_name, min_q),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue_d)[1]/2,height = 5+length(wanted_labels)/10)
#   
# }

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

## setting new metabololmics data
map2 = read.table(file = '../../metabolomics/China/data/SOURCE_china_mets_table_norm_v2.tsv', 
                  header = T, row.names = 1, sep = '\t', quote = '',comment.char = '')
map2 = as.data.frame(t(map2))
metadata_df2 = data.frame(SampleID = row.names(metadata_df), 
                          # pn_ID = gsub('TI.*','',row.names(metadata_df)))
                          pn_ID = gsub(sprintf('%s.*',tt),'',row.names(metadata_df)))
map2$pn_ID = row.names(map2)
metadata_df2 = merge(metadata_df2, map2, by = 'pn_ID',all.x = T)
row.names(metadata_df2) = metadata_df2$SampleID
metadata_df2 = metadata_df2[,c(-1,-2)]

# log scaleong and cleaning metaboloimcs data (same as in Israel data)
met = metadata_df2
added_data_type = 'LogClean'
if( grepl('LogClean',added_data_type) )
{
  for ( i in 1:dim(met)[2] )
  {
    # change zeroes to 5th of the smallest value per metabolite
    met[met[,i] == 0 & !is.na(met[,i]),i] = min(met[met[,i] != 0& !is.na(met[,i]),i], na.rm = T)/5
    # change values over 4 SDs from the mean to the top/bottom value without it
    top_cut = mean(met[,i], na.rm = T) + 4*sd(met[,i], na.rm = T)
    met[met[,i] > top_cut& !is.na(met[,i]),i ] = max(met[ met[,i]<=top_cut & !is.na(met[,i]),i ], na.rm = T)
    bottom_cut = mean(met[,i], na.rm = T) - 4*sd(met[,i], na.rm = T)
    met[met[,i] < bottom_cut & !is.na(met[,i]),i ] = min(met[ met[,i]>=bottom_cut & !is.na(met[,i]),i ], na.rm = T)
  }
  # log transform
  met = log10(met)
}
metadata_df2 = met

## ME calculation
# Calculate MEs with color labels
# name = sprintf('SOURCE_China_TI_by_%s', israel_name)
name = sprintf('SOURCE_China_ageMatch_TI_by_%s', israel_name)
# WGCNA_plot_genes_correlaitons_by_module(ftr_df, israel_df$module_color, out_path)
res = WGCNA_calculate_module_eigengenes_and_corrs(ftr_df,israel_df$module_color, metadata_df2, cor_type, name)
MEs = res[[1]]
moduleTraitCor = res[[2]]
moduleTraitPvalue = res[[3]]


# my version
heatmap = WGCNA_ME_met_heatmap(moduleTraitCor, moduleTraitPvalue, name, metadata_prms)
ggsave(sprintf('%s/%s_module_cor_heatmap.pdf',out_path, name),plot = heatmap, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+dim(moduleTraitPvalue)[2]/10)
# if (plot_disease_related_modules_flag)
# {
#   heatmap2 = WGCNA_ME_met_heatmap_prm_related(moduleTraitCor, moduleTraitPvalue, name, prm = 'Disease', cutoff = min_p, metadata_prms)
#   ggsave(sprintf('%s/%s_module_cor_heatmap_%sF%s.pdf',out_path, name,'Disease',min_p), plot = heatmap2, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+dim(moduleTraitPvalue)[2]/10)
# }


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
mdc2 = res_china[[2]]

heatmap_s_all = ggplot(mdc2, aes(x=module,y=variable,fill = correlation)) +
  geom_tile() + xlab('Module') +
  scale_fill_gradientn(colours = col2, limits = c(-1,1)) +
  # geom_point(aes(shape = p_value<=0.05)) + scale_shape_manual(values = c(NA, 8)) +
  # geom_text(aes(label=sprintf('%.2f\n%.2e',correlation,p_value )), size=3) +
  # geom_text(aes(label=sprintf('%.2e',p_value )), size=3) +
  scale_colour_manual(values = c('black','white')) +
  geom_text(aes(label=signif(p_value, 1), colour =p_value < 0.001 ), size=3) +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size= 10)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size= 10)) + coord_flip() +
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  ggtitle('China metabolomics all') +
  # theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  ylab('Metabolite')
y_lables = ggplot_build(heatmap_s)$layout$panel_params[[1]]$y$get_labels()
ggsave(sprintf('%s/%s_module_cor_heatmap_q%s.pdf',out_path2, name, min_q),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+length(y_lables)/10)

write.table(x = mdc2, file = sprintf('%s/%s_module_cor_df.txt', out_path2, name), quote = F, sep = '\t', row.names = F) 

# heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor[wanted_modules,], moduleTraitPvalue[wanted_modules,], name, metadata_prms, min_p = min_p, wanted_modules_order = wanted_modules, fdr_flag_non_metadata = F,  show_p_value_fdr_flag = T)
# y_lables = ggplot_build(heatmap_s)$layout$panel_params[[1]]$y$get_labels()
# ggsave(sprintf('%s/%s_module_cor_heatmap_p%s.pdf',out_path2, name, min_p),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+length(y_lables)/10)
# 
if (subcohorts_flag)
{
  # for looking at only a subcohort ( only diease/healthy)
  subcohort_name = 'Disease'
  wanted_pos = metadata_df$Disease == 1

  res = WGCNA_calculate_module_eigengenes_and_corrs(ftr_df[wanted_pos,],israel_df$module_color, clean_metadata_eampty_noVar( metadata_df2[wanted_pos,] ), cor_type, name)
  moduleTraitCor_d = res[[2]]
  moduleTraitPvalue_d = res[[3]]
  wanted_labels = y_lables[y_lables %in% colnames(moduleTraitCor_d)]

  heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor_d[wanted_modules,], moduleTraitPvalue_d[wanted_modules,], name, metadata_prms, wanted_modules_order = wanted_modules)
  # ggsave(sprintf('%s/%s_module_cor_heatmap_%s.pdf',out_path2, name,subcohort_name),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue_d)[1]/2,height = 5+dim(moduleTraitPvalue_d)[2]/10)

  res = WGCNA_ME_met_heatmap_full_res(moduleTraitCor_d[wanted_modules,wanted_labels], moduleTraitPvalue_d[wanted_modules,wanted_labels], name, metadata_prms, min_p = 1, wanted_modules_order = wanted_modules, fdr_flag_non_metadata=  T, cluster_cols_order_flag = F, show_p_value_fdr_flag = F)
  heatmap_s_cd = res[[1]]
  mdc2_cd = res[[2]]
  heatmap_s_cd = WGCNA_ME_met_heatmap(moduleTraitCor_d[wanted_modules,wanted_labels], moduleTraitPvalue_d[wanted_modules,wanted_labels], name, metadata_prms, wanted_modules_order = wanted_modules, cluster_cols_order_flag = F )
  ggsave(sprintf('%s/%s_module_cor_heatmap_%s_allq%s.pdf',out_path2, name,subcohort_name, min_q),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue_d)[1]/2,height = 5+length(wanted_labels)/10)

  subcohort_name = 'Healthy'
  wanted_pos = metadata_df$Disease == 0

  res = WGCNA_calculate_module_eigengenes_and_corrs(ftr_df[wanted_pos,],israel_df$module_color, clean_metadata_eampty_noVar( metadata_df2[wanted_pos,] ), cor_type, name)
  moduleTraitCor_d = res[[2]]
  moduleTraitPvalue_d = res[[3]]
  wanted_labels = y_lables[y_lables %in% colnames(moduleTraitCor_d)]

  heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor_d[wanted_modules,], moduleTraitPvalue_d[wanted_modules,], name, metadata_prms, wanted_modules_order = wanted_modules)
  # ggsave(sprintf('%s/%s_module_cor_heatmap_%s.pdf',out_path2, name,subcohort_name),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue_d)[1]/2,height = 5+dim(moduleTraitPvalue_d)[2]/10)

  heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor_d[wanted_modules,wanted_labels], moduleTraitPvalue_d[wanted_modules,wanted_labels], name, metadata_prms, wanted_modules_order = wanted_modules, cluster_cols_order_flag = F )
  ggsave(sprintf('%s/%s_module_cor_heatmap_%s_allq%s.pdf',out_path2, name,subcohort_name, min_q),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue_d)[1]/2,height = 5+length(wanted_labels)/10)


  # g = ggplot(metadata_df, aes(x=as.factor(metadata_df$Active_Smoker_manual), y=MEs$MEtan)) +
  #   geom_boxplot(outlier.alpha = 0) + geom_jitter(aes(colour = as.factor(metadata_df$Disease))) +
  #   scale_color_manual(values = c('blue','red'), name = 'Disease') + xlab('Active_Smoker_manual') + # ylab('tan') +
  #   theme_bw()
  #
  # ggsave(sprintf('%s/%s_smoking_active_tan_box.pdf',out_path2, name),plot = g, device = 'pdf', width = 4,height = 3.5)
}
p  = ( heatmap_s_all + theme(legend.position = 'none') + 
         ggtitle('China all CD and controls') ) + 
          ( heatmap_s_cd + theme( axis.text.y = element_blank(), axis.ticks.y = element_blank() ) + 
         ggtitle('China only within CD') + guides(color = "none") )
out_path = 'rnaSeq_main/final_figs/'
ggsave(sprintf('%s/module_cor_heatmap_china_mets_all_and_CD_allq%s.pdf',out_path, min_q),plot = p, device = 'pdf', width = 12,height = 25)

mdc2$group = 'CD and control'
mdc2_cd$group = 'CD'
mdc2$Country = 'China'
mdc2$var_type = 'Metabolites'
mdc2_cd$Country = 'China'
mdc2_cd$var_type = 'Metabolites'
temp = rbind(mdc2, mdc2_cd)
write.table(x = temp,
            file = sprintf('%s/module_cor_heatmap_china_mets_q%s.txt',out_path , min_q),
            quote = F, sep = '\t', row.names = F)
mdc = mdc2
mdc2 = mdc2_cd
mdc$id = sprintf('%s_%s',mdc$module, mdc$variable)
mdc2$id = sprintf('%s_%s',mdc2$module, mdc2$variable)
mdc = mdc[order(mdc$id),]
mdc2 = mdc2[order(mdc2$id),]

t = mdc
t$cd_cor = mdc2$correlation
t=t[t$q_value<=0.25,]
t$same_dir = ifelse( (t$correlation < 0 & t$cd_cor < 0) | (t$correlation > 0 & t$cd_cor > 0),T,F )
print( sum(t$same_dir)/length(t$same_dir) )
print(table(t$same_dir))

print(table(t$cd_cor>0))
p_cd = sum( t$cd_cor>0) / length(t$cd_cor) 
print(table(t$correlation>0))
p_all = sum( t$correlation>0) / length(t$correlation) 
prob = p_cd*p_cd + (1-p_cd)*(1-p_all)
binom.test(x = sum(t$same_dir==T), n = length(t$cd_cor) , p = prob)




