library(ggplot2)

## filter metabolomics data using WGCNA modules
modules_mets_file = '../WGCNA/rnaSeq_main/res_v2/SOURCE_Israel_TI_withStoolMetablolomicsV2LogClean_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default//selected_modules//SOURCE_Israel_TI_withStoolMetablolomicsV2LogClean_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default_module_cor_df.txt'
md_mets = read.table(file = modules_mets_file, header = T, sep = '\t')
md_mets$variable = make.names(md_mets$variable)
metadata_prms = c('Gender','Age_years','BMI','Disease','inflammed','CRP_numeric','Calprotectin_numeric','Active_Smoker_manual','Previous_or_Current_Tobacco_Use_manual','health_index_same_sample','health_index_stool')
md_mets = md_mets[! md_mets$variable %in% c(metadata_prms),]

wanted_modules = c('yellow','green','red','pink', 'purple','tan','salmon','black','brown')
res_df = data.frame(module = wanted_modules)
res_df_dirs = data.frame(module = wanted_modules)

dx_down_modules = c('yellow','green','red','pink')
dx_up_modules = c('purple','tan','salmon','black','brown')

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
  
  FFQ_halla_file = sprintf('res/full_pairs/Israel_metabolomicsM%s_Stool_vs_Israel_FFQ/halla_res_alpha0.25/all_associations.txt',md)
  FFQ_halla = read.table(file = FFQ_halla_file, header = T, sep = '\t')
  FFQ_halla = FFQ_halla[ FFQ_halla$q.values <= 0.25, ]
  res_df$halla_ffq_n[i] = length( FFQ_halla$X_features )
  res_df_dirs$halla_ffq_n_pos[i] = length( FFQ_halla$X_features[FFQ_halla$association > 0] )
  res_df_dirs$halla_ffq_n_neg[i] = length( FFQ_halla$X_features[FFQ_halla$association < 0] )
  
  taxa_halla_file = sprintf('res/full_pairs/Israel_metabolomicsM%s_Stool_vs_Israel_16S_Stool/halla_res_alpha0.25/all_associations.txt',md)
  taxa_halla = read.table(file = taxa_halla_file, header = T, sep = '\t')
  taxa_halla = taxa_halla[ taxa_halla$q.values <= 0.25, ]
  res_df$halla_16S_n[i] = length( taxa_halla$X_features )
  res_df_dirs$halla_16S_n_pos[i] = length( taxa_halla$X_features[taxa_halla$association > 0] )
  res_df_dirs$halla_16S_n_neg[i] = length( taxa_halla$X_features[taxa_halla$association < 0] )
  
  taxa_halla_file = sprintf('res/full_pairs/Israel_metabolomicsM%s_Stool_vs_Israel_metagenomics_pathabundance/halla_res_alpha0.25/all_associations.txt',md)
  taxa_halla = read.table(file = taxa_halla_file, header = T, sep = '\t')
  taxa_halla = taxa_halla[ taxa_halla$q.values <= 0.25, ]
  res_df$halla_metagenomics_pathabundance_n[i] = length( taxa_halla$X_features )
  res_df_dirs$halla_metagenomics_pathabundance_n_pos[i] = length( taxa_halla$X_features[taxa_halla$association > 0] )
  res_df_dirs$halla_metagenomics_pathabundance_n_neg[i] = length( taxa_halla$X_features[taxa_halla$association < 0] )
  
  taxa_halla_file = sprintf('res/full_pairs/Israel_metabolomicsM%s_Stool_vs_Israel_metagenomics_ecs/halla_res_alpha0.25/all_associations.txt',md)
  taxa_halla = read.table(file = taxa_halla_file, header = T, sep = '\t')
  taxa_halla = taxa_halla[ taxa_halla$q.values <= 0.25, ]
  res_df$halla_metagenomics_ecs_n[i] = length( taxa_halla$X_features )
  res_df_dirs$halla_metagenomics_ecs_n_pos[i] = length( taxa_halla$X_features[taxa_halla$association > 0] )
  res_df_dirs$halla_metagenomics_ecs_n_neg[i] = length( taxa_halla$X_features[taxa_halla$association < 0] )
  
  taxa_halla_file = sprintf('res/full_pairs/Israel_metabolomicsM%s_Stool_vs_Israel_metagenomics_species/halla_res_alpha0.25/all_associations.txt',md)
  taxa_halla = read.table(file = taxa_halla_file, header = T, sep = '\t')
  taxa_halla = taxa_halla[ taxa_halla$q.values <= 0.25, ]
  res_df$halla_metagenomics_species_n[i] = length( taxa_halla$X_features )
  res_df_dirs$halla_metagenomics_species_n_pos[i] = length( taxa_halla$X_features[taxa_halla$association > 0] )
  res_df_dirs$halla_metagenomics_species_n_neg[i] = length( taxa_halla$X_features[taxa_halla$association < 0] )
  
}

res_df_m = reshape2::melt(res_df)
res_df_dirs_m = reshape2::melt(res_df_dirs)

res_df_m$module = factor(res_df_m$module, levels = wanted_modules)
res_df_m$module_direction = ifelse(res_df_m$module %in% dx_up_modules, 'With disease','Against disease')
res_df_m$variable2 = gsub('halla_','',res_df_m$variable)
res_df_m$variable2 = gsub('_n','',res_df_m$variable2)
res_df_m$variable2 = gsub('_','\n',res_df_m$variable2)

res_df_mf = res_df_m[grepl(pattern = 'n_up|n_down', res_df_m$variable),]
res_df_mf$module_direction = ifelse(res_df_mf$module %in% dx_up_modules, 'Up in disease','Up in control')
res_df_mf$module_direction = factor(res_df_mf$module_direction, levels = unique(res_df_mf$module_direction))
res_df_mf$variable2 = gsub('_up_|_down_','_',res_df_mf$variable)
res_df_mf$value[res_df_mf$value==0]=NA
g2_sep2_n = ggplot(res_df_mf, 
                   aes(x=module, y=value, group = variable, fill = variable)) + 
  geom_bar(stat="identity", position=position_dodge(), colour = 'black', width = 0.8) + 
  theme_classic() +
  # scale_fill_brewer(name = '', palette = 'Blues' ) + 
  scale_fill_manual(name = '', values = c('#FFA372','#7FDBDA','#4E89AE','#ED6663'), 
                    breaks = c('module_corr_mets_n_up_pos','module_corr_mets_n_up_neg',
                               'module_corr_mets_n_down_pos','module_corr_mets_n_down_neg'),  
                    labels = c('Disease positive','Disease negative','Control positive','Control negative' ) 
  ) + 
  ylab('N of FDR < 0.25\nfor metabolites correaltions\nwith WGCNA modules') +
  xlab('Module') + # facet_grid(variable~.) +
  facet_grid(.~module_direction, scales = 'free_x') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        # strip.background = element_blank(),
        # strip.text = element_blank()
  ) 
ggsave('res/n_WGCNA_modules_realted_mets2.pdf',plot = g2_sep2_n, device = 'pdf', width = 4.5,height = 2.5, limitsize = FALSE, scale = 1.2)


temp = res_df_dirs_m[grepl(pattern = '_pos|_neg', res_df_dirs_m$variable),]
temp$module_direction = ifelse(temp$module %in% dx_up_modules, 'Up in disease','Up in control')
temp$module_direction = factor(temp$module_direction, levels = unique(temp$module_direction))
temp$variable2 = gsub('_pos|_neg','',temp$variable)
temp$variable2 = gsub('halla_','',temp$variable2)
temp$variable2 = gsub('_n','',temp$variable2)
temp$variable2 = gsub('_','\n',temp$variable2)
temp$dir = gsub('.*_','',temp$variable)
temp$value[temp$value==0]=NA
temp$module = factor(temp$module, levels = wanted_modules)
temp$variable2 = gsub('ffq','FFQ',temp$variable2)
temp$variable2 = factor(temp$variable2, levels = c('FFQ','16S',
                                                   'metagenomics\nspecies',
                                                   'metagenomics\npathabundance',
                                                   'metagenomics\necs'))
g3_3 = ggplot(temp, 
              aes(x=module, y=value, group = dir, fill = dir)) + 
  # aes(x=module, y=value, group = dir, fill = variable2, alpha = dir)) + 
  geom_bar(stat="identity", position=position_stack(), colour = 'black', width = 0.8) + 
  theme_classic() +
  # guides(fill = 'none') + 
  # scale_fill_manual(values = c('skyblue','#EB455F'), name= 'Correlation\ndireciton') + 
  # scale_fill_manual(values = c('gray75','gray45'), name= 'Correlation\ndireciton') + 
  # scale_fill_manual(values = c('gray60','gray60', 'gray60'), ) + 
  scale_fill_manual(values = c('#ebf1f5','#769cbc'),name= 'Correlation\ndireciton' ) + 
  # scale_fill_brewer(name = '', palette = 'Blues' ) + 
  # scale_alpha_manual(values = c(0.3,1), name= 'Correlation\ndireciton') + 
  # scale_fill_manual(name = '', values = c('#FFA372','#7FDBDA','#4E89AE','#ED6663'), 
  #                   breaks = c('module_corr_mets_n_up_pos','module_corr_mets_n_up_neg',
  #                              'module_corr_mets_n_down_pos','module_corr_mets_n_down_neg'),  
  #                   labels = c('Disease positive','Disease negative','Control positive','Control negative' ) 
  # ) + 
  # ylab('N of FDR < 0.25\nfor metabolites correaltions\nwith WGCNA modules') +
  ylab('Number of significant\nassociations') +
  xlab('Module') + # facet_grid(variable~.) +
  facet_grid(variable2~module_direction, scales = 'free') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        strip.text.y = element_text(angle = 0),
        # legend.position ='none'
        # strip.background = element_blank(),
        # strip.text = element_blank()
  ) 
ggsave('res/n_corr_WGCNA_modules_realted_mets_pos_neg_3.pdf',plot = g3_3, device = 'pdf', 
       width = 6,height = 4, limitsize = FALSE)



