library(ggplot2)

## filter metabolomics data using WGCNA modules
# modules_mets_file = '../WGCNA/rnaSeq_main/res_v2/SOURCE_Israel_TI_withStoolMetablolomicsV2LogClean_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default//selected_modules//SOURCE_Israel_TI_withStoolMetablolomicsV2LogClean_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default_module_cor_df.txt'
modules_mets_file = '../WGCNA/rnaSeq_main/res_other_modules/TI_china_israel_modules_metsLogClean/selected_modules/SOURCE_China_ageMatch_TI_by_SOURCE_Israel_TI_withFFQv10_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default_module_cor_df.txt'
md_mets = read.table(file = modules_mets_file, header = T, sep = '\t')
md_mets$variable = make.names(md_mets$variable)
metadata_prms = c('Gender','Age_years','BMI','Disease','inflammed','CRP_numeric','Calprotectin_numeric','Active_Smoker_manual','Previous_or_Current_Tobacco_Use_manual','health_index_same_sample','health_index_stool')
md_mets = md_mets[! md_mets$variable %in% c(metadata_prms),]

wanted_modules = c('yellow','green','red','pink', 'purple','tan','salmon','black','brown')
res_df = data.frame(module = wanted_modules)
res_df_dirs = data.frame(module = wanted_modules)

dx_down_modules = c('yellow','green','red','pink')
dx_up_modules = c('purple','tan','salmon','black','brown')

# halla_path = 'res/full_pairs/'
halla_path = 'res/full_pairs_china_mds_mets/'
group_add = ''

# halla_path = 'res/full_pairs_china_mds_mets_wgcna_samples/'
# group_add = '_WGCNAsamps'

# halla_path = 'res/full_pairs_china_mds_mets_groups/'
# # group_add = '_rural'
# group_add = '_urban_CD'
# # group_add = '_urban'
# # group_add = '_CD'

for ( i in 1:length(wanted_modules) )
{
  md = wanted_modules[i]
  
  up_mets = md_mets$variable[ md_mets$q_value <= 0.25 & md_mets$module == md & md_mets$correlation > 0]
  down_mets = md_mets$variable[ md_mets$q_value <= 0.25 & md_mets$module == md & md_mets$correlation < 0]
  mets = md_mets$variable[ md_mets$q_value <= 0.25 & md_mets$module == md]
  res_df$module_corr_mets_n[i] = length(mets)
  res_df$module_corr_mets_n_pos[i] = length(up_mets)
  res_df$module_corr_mets_n_neg[i] = length(down_mets)
  
  up_pos_mets = md_mets$variable[ md_mets$q_value <= 0.25 & md_mets$module == md & md_mets$correlation > 0 & md %in% dx_up_modules]
  up_neg_mets = md_mets$variable[ md_mets$q_value <= 0.25 & md_mets$module == md & md_mets$correlation < 0 & md %in% dx_up_modules]
  down_pos_mets = md_mets$variable[ md_mets$q_value <= 0.25 & md_mets$module == md & md_mets$correlation > 0 & md %in% dx_down_modules]
  down_neg_mets = md_mets$variable[ md_mets$q_value <= 0.25 & md_mets$module == md & md_mets$correlation < 0 & md %in% dx_down_modules]
  res_df$module_corr_mets_n_up_pos[i] = length(up_pos_mets)
  res_df$module_corr_mets_n_up_neg[i] = length(up_neg_mets)
  res_df$module_corr_mets_n_down_pos[i] = length(down_pos_mets)
  res_df$module_corr_mets_n_down_neg[i] = length(down_neg_mets)
  
  FFQ_halla_file = sprintf('%s/China_metabolomicsM%s_Stool_vs_China_FFQ%s/halla_res_alpha0.25/all_associations.txt',halla_path, md, group_add)
  FFQ_halla = read.table(file = FFQ_halla_file, header = T, sep = '\t')
  FFQ_halla = FFQ_halla[ FFQ_halla$q.values <= 0.25, ]
  res_df$halla_ffq_n[i] = length( FFQ_halla$X_features )
  res_df_dirs$halla_ffq_n_pos[i] = length( FFQ_halla$X_features[FFQ_halla$association > 0] )
  res_df_dirs$halla_ffq_n_neg[i] = length( FFQ_halla$X_features[FFQ_halla$association < 0] )
  
  taxa_halla_file = sprintf('%s/China_metabolomicsM%s_Stool_vs_China_16S_Stool%s/halla_res_alpha0.25/all_associations.txt',halla_path ,md, group_add)
  taxa_halla = read.table(file = taxa_halla_file, header = T, sep = '\t')
  taxa_halla = taxa_halla[ taxa_halla$q.values <= 0.25, ]
  res_df$halla_16S_n[i] = length( taxa_halla$X_features )
  res_df_dirs$halla_16S_n_pos[i] = length( taxa_halla$X_features[taxa_halla$association > 0] )
  res_df_dirs$halla_16S_n_neg[i] = length( taxa_halla$X_features[taxa_halla$association < 0] )
  
  taxa_halla_file = sprintf('%s/China_metabolomicsM%s_Stool_vs_China_metagenomics_pathabundance%s/halla_res_alpha0.25/all_associations.txt',halla_path,md, group_add)
  taxa_halla = read.table(file = taxa_halla_file, header = T, sep = '\t')
  taxa_halla = taxa_halla[ taxa_halla$q.values <= 0.25, ]
  res_df$halla_metagenomics_pathabundance_n[i] = length( taxa_halla$X_features )
  res_df_dirs$halla_metagenomics_pathabundance_n_pos[i] = length( taxa_halla$X_features[taxa_halla$association > 0] )
  res_df_dirs$halla_metagenomics_pathabundance_n_neg[i] = length( taxa_halla$X_features[taxa_halla$association < 0] )

  taxa_halla_file = sprintf('%s/China_metabolomicsM%s_Stool_vs_China_metagenomics_ecsNamed%s/halla_res_alpha0.25/all_associations.txt',halla_path,md, group_add)
  taxa_halla = read.table(file = taxa_halla_file, header = T, sep = '\t')
  taxa_halla = taxa_halla[ taxa_halla$q.values <= 0.25, ]
  res_df$halla_metagenomics_ecs_n[i] = length( taxa_halla$X_features )
  res_df_dirs$halla_metagenomics_ecs_n_pos[i] = length( taxa_halla$X_features[taxa_halla$association > 0] )
  res_df_dirs$halla_metagenomics_ecs_n_neg[i] = length( taxa_halla$X_features[taxa_halla$association < 0] )

  taxa_halla_file = sprintf('%s/China_metabolomicsM%s_Stool_vs_China_metagenomics_species%s/halla_res_alpha0.25/all_associations.txt',halla_path,md, group_add)
  taxa_halla = read.table(file = taxa_halla_file, header = T, sep = '\t')
  taxa_halla = taxa_halla[ taxa_halla$q.values <= 0.25, ]
  res_df$halla_metagenomics_species_n[i] = length( taxa_halla$X_features )
  res_df_dirs$halla_metagenomics_species_n_pos[i] = length( taxa_halla$X_features[taxa_halla$association > 0] )
  res_df_dirs$halla_metagenomics_species_n_neg[i] = length( taxa_halla$X_features[taxa_halla$association < 0] )

}

res_df_m = reshape2::melt(res_df)
res_df_dirs_m = reshape2::melt(res_df_dirs)

res_df_m$module = factor(res_df_m$module, levels = wanted_modules)
g = ggplot(res_df_m[res_df_m$variable =='module_corr_mets_n',], aes(x=module, y=value, group = variable, fill = variable)) + 
  geom_bar(stat="identity", position=position_dodge(), colour = 'black', width = 0.8) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab('N of FDR < 0.25\nfor metabolites correlated with\nrnaSeq TI modules') + theme(legend.position = 'none') + 
  xlab('Module') + 
  scale_fill_brewer(name = '', palette = 'Paired' ) 
# ggsave('res/n_metabolites_corr_WGCNA_modules.pdf',plot = g, device = 'pdf', width = 5,height = 3, limitsize = FALSE)
g_up_down = ggplot(res_df_m[res_df_m$variable %in% c('module_corr_mets_n_pos','module_corr_mets_n_neg'),], aes(x=module, y=value, group = variable, fill = variable)) + 
  geom_bar(stat="identity", position=position_dodge(), colour = 'black', width = 0.8) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab('N of FDR < 0.25\nfor metabolites correlated with\nrnaSeq TI modules') + 
  # theme(legend.position = 'none') + 
  xlab('Module') + 
  scale_fill_brewer(name = '', palette = 'Paired' ) 
# ggsave('res/n_metabolites_corr_WGCNA_modules_cor_direction.pdf',plot = g_up_down, device = 'pdf', width = 6,height = 3, limitsize = FALSE)

# g2 = ggplot(res_df_m[res_df_m$variable !='module_corr_mets_n' & !grepl('metagenomic',res_df_m$variable),], aes(x=module, y=value, group = variable, fill = variable)) + 
# g2 = ggplot(res_df_m[res_df_m$variable %in% c('halla_16S_n','halla_ffq_n','halla_metagenomics_species_n'),], aes(x=module, y=value, group = variable, fill = variable)) + 
g2 = ggplot(res_df_m[res_df_m$variable %in% c('halla_16S_n','halla_ffq_n'),], aes(x=module, y=value, group = variable, fill = variable)) + 
  geom_bar(stat="identity", position=position_dodge(), colour = 'black', width = 0.8) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab('N of FDR < 0.25\nfor halla correaltions\nwith metabolites') +
  xlab('Module') # + 
# scale_fill_brewer(name = '', labels = c('FFQ','16S'),palette = 'Paired' ) 
# scale_fill_brewer(name = '', palette = 'Paired' )# +scale_y_log10()
# ggsave('res/n_FFQ_16S_corr_WGCNA_modules_realted_mets.pdf',plot = g2, device = 'pdf', width = 6,height = 3, limitsize = FALSE)

res_df_m$module_direction = ifelse(res_df_m$module %in% dx_up_modules, 'With disease','Against disease')
res_df_m$variable2 = gsub('halla_','',res_df_m$variable)
res_df_m$variable2 = gsub('_n','',res_df_m$variable2)
res_df_m$variable2 = gsub('_','\n',res_df_m$variable2)
g2_sep_corrs_temp  = res_df_m[! grepl(pattern = 'mets_n', res_df_m$variable),]
# g2_sep_corrs_temp$variable2 = factor(g2_sep_corrs_temp$variable2, levels = c('ffq','16S'))
g2_sep_corrs_temp$variable2 = factor(g2_sep_corrs_temp$variable2, levels = c('ffq','16S',
                                                                             'metagenomics\nspecies',
                                                                             'metagenomics\necs',
                                                                             'metagenomics\npathabundance'))
g2_sep_corrs = ggplot(g2_sep_corrs_temp, 
                      aes(x=module, y=value, group = variable2, fill = variable2)) + 
  geom_bar(stat="identity", position=position_dodge(), colour = 'black', width = 0.6) + 
  theme_classic() +
  scale_fill_brewer(name = '', palette = 'Blues' ) + 
  ylab('Number of significant associations') +
  xlab('Module') + 
  facet_grid(variable2~module_direction, scales = 'free') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        # strip.background.y = element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position ='none'
  ) 
ggsave(sprintf('res/china_n_corr_WGCNA_modules_realted_mets%s.pdf',group_add ),plot = g2_sep_corrs, device = 'pdf', width = 5,height = 5, limitsize = FALSE)

temp = res_df_m[! grepl(pattern = 'mets_n', res_df_m$variable) &
                  ! grepl(pattern = 'ecs', res_df_m$variable) & 
                  ! grepl(pattern = 'pathabundance', res_df_m$variable),]
temp$variable2 = factor(temp$variable2, levels = c('ffq','16S','metagenomics\nspecies'))
temp$variable3 = as.character(temp$variable2)
temp$variable3[temp$variable3=='ffq'] = 'FFQ'
temp$variable3[temp$variable3=='16S'] = '16S\ntaxonomy'
temp$variable3[temp$variable3=='metagenomics\nspecies'] = 'MGX\nspecies'
temp$variable3 = factor(temp$variable3, levels = c('FFQ','16S\ntaxonomy','MGX\nspecies'))

g2_sep_corrs_f = ggplot(temp, 
                        aes(x=module, y=value, group = variable2, fill = variable2)) + 
  geom_bar(stat="identity", position=position_dodge(), colour = 'black', width = 0.6) + 
  theme_classic() +
  scale_fill_brewer(name = '', palette = 'Blues' ) + 
  # ylab('N of FDR < 0.25\nfor halla correaltions\nwith metabolites') +
  ylab('Number of significant associations') +
  xlab('Module') + facet_grid(variable3~module_direction, scales = 'free') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        # strip.background.y = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0),
        legend.position ='none',strip.background=element_rect(color = NA)
  ) 
# ggsave('res/n_corr_WGCNA_modules_realted_mets_f.pdf',plot = g2_sep_corrs_f, device = 'pdf', width = 5,height = 3, limitsize = FALSE)

sum(temp$value[temp$variable3 == 'FFQ' & temp$module_direction == 'With disease'])
sum(temp$value[temp$variable3 == 'FFQ' & temp$module_direction == 'Against disease'])
sum(temp$value[temp$variable3 == '16S\ntaxonomy' & temp$module_direction == 'With disease'])
sum(temp$value[temp$variable3 == '16S\ntaxonomy' & temp$module_direction == 'Against disease'])
sum(temp$value[temp$variable3 == 'MGX\nspecies' & temp$module_direction == 'With disease'])
sum(temp$value[temp$variable3 == 'MGX\nspecies' & temp$module_direction == 'Against disease'])

