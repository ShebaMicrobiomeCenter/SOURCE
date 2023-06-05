library(WGCNA)
library(ggplot2)
options(stringsAsFactors = FALSE);

source('/pita/users/tzipi/code/R_figs/WGCNA_funcs.R')

cor_type = 'WGCNA_default'

in_path = 'rnaSeq_main/res_v2/'
out_path = 'rnaSeq_main/res_other_modules/'
out_path = sprintf('%s/TI_israel_pcm/', out_path)
israel_name = sprintf('SOURCE_Israel_TI_withFFQv10_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_%s',cor_type)
# china_name = sprintf('SOURCE_China_TI_ageMatch_withFFQv10_oneSample_pwr_14_netType_signed_hybrid_minMdSz_30_%s',cor_type)
dir.create(out_path)

wanted_modules = c('yellow','green','red','pink', 'purple','tan','salmon','black','brown')
wanted_modules = sprintf('ME%s',wanted_modules)
out_path2 = sprintf('%s/selected_modules/', out_path)
dir.create(out_path2)

name = israel_name

israel_ti = load( sprintf('%s/%s/data.RData',in_path, name) )

metadata_df 

metadata_prms = c('Gender','Age_years','BMI',
                  'Active_Smoker_manual','Previous_or_Current_Tobacco_Use_manual',
                  'Disease','inflammed','CRP_numeric','Calprotectin_numeric',# 'CRP_log','Calprotectin_log',
                  'health_index_same_sample','health_index_stool')

metadata_df = metadata_df[,names(metadata_df) %in%metadata_prms ]
new_data_file = '/pita/users/tzipi/projects/multiomics/SOURCE/metadata/sheba_PCM.txt'
ndf = read.table(new_data_file, header = T, sep = '\t')
ndf$type = ifelse(ndf$type == 'Type_I',1,2)
metadata_ids = row.names(metadata_df)
metadata_ids = gsub('TI[12]|R[34]','',metadata_ids)
metadata_ids = gsub('_dup1$','',metadata_ids)
metadata_ids = gsub('[.]$','',metadata_ids)
metadata_ids = gsub('.*[.]','',metadata_ids)
metadata_df$pn_ID = metadata_ids
metadata_df2 = merge(metadata_df, ndf, by='pn_ID',all.x = T)
row.names(metadata_df2) = metadata_df2$pn_ID
metadata_df2 = metadata_df2[metadata_df$pn_ID,]
row.names(metadata_df2) = row.names(metadata_df)
metadata_df = metadata_df2[,-1]

res = WGCNA_calculate_module_eigengenes_and_corrs(ftr_df,df$module_color, metadata_df, cor_type, name)
MEs = res[[1]]
moduleTraitCor = res[[2]]
moduleTraitPvalue = res[[3]]
heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor[wanted_modules,], moduleTraitPvalue[wanted_modules,], name, metadata_prms, wanted_modules_order = wanted_modules)
# ggsave(sprintf('%s/%s_module_cor_heatmap_pcm.pdf',out_path2, name),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+dim(moduleTraitPvalue)[2]/10)

subcohort_name = 'Disease'
wanted_pos = metadata_df$Disease == 1

res = WGCNA_calculate_module_eigengenes_and_corrs(ftr_df[wanted_pos,],df$module_color, clean_metadata_eampty_noVar( metadata_df[wanted_pos,] ), cor_type, name)
MEs = res[[1]]
moduleTraitCor = res[[2]]
moduleTraitPvalue = res[[3]]
heatmap_s = WGCNA_ME_met_heatmap(moduleTraitCor[wanted_modules,], moduleTraitPvalue[wanted_modules,], name, metadata_prms, wanted_modules_order = wanted_modules)
# ggsave(sprintf('%s/%s_module_cor_heatmap_pcm_disease.pdf',out_path2, name),plot = heatmap_s, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+dim(moduleTraitPvalue)[2]/10)

