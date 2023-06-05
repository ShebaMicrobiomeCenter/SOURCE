source('/pita/users/tzipi/code/R_figs/figs_funcs.R')

# out_path  = 'rnaSeq_main/module_related_mets/'
out_path  = 'rnaSeq_main/module_related_mets_china_ti/'
dir.create(out_path)

## filter metabolomics data using WGCNA modules
modules_mets_file = 'rnaSeq_main/res_v2/SOURCE_Israel_TI_withStoolMetablolomicsV2LogClean_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default//selected_modules//SOURCE_Israel_TI_withStoolMetablolomicsV2LogClean_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default_module_cor_df.txt'
modules_mets_file = 'rnaSeq_main/res_other_modules/TI_china_israel_modules_metsLogClean/selected_modules/SOURCE_China_ageMatch_TI_by_SOURCE_Israel_TI_withFFQv10_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default_module_cor_df.txt'
md_mets = read.table(file = modules_mets_file, header = T, sep = '\t')
md_mets$variable = make.names(md_mets$variable)
metadata_prms = c('Gender','Age_years','BMI','Disease','inflammed','CRP_numeric','Calprotectin_numeric','Active_Smoker_manual','Previous_or_Current_Tobacco_Use_manual','health_index_same_sample','health_index_stool')
md_mets = md_mets[! md_mets$variable %in% c(metadata_prms),]

md_mets_all = md_mets
md_mets = md_mets[md_mets$q_value <= 0.25,]

wanted_modules = c('yellow','green','red','pink', 'purple','tan','salmon','black','brown')
for (md in wanted_modules)
{
  md_mets = md_mets_all
  # 
  md_mets = md_mets[md_mets$q_value <= 0.25 & md_mets$module == md,]
  # name1 = gsub(pattern = 'metabolomics',sprintf('metabolomicsM%s',md),name1)
  ## 
  up = md_mets$variable[md_mets$correlation > 0]
  down = md_mets$variable[md_mets$correlation < 0]
  
  write.table(x = up, file = sprintf('%s/%s_up.txt',out_path, md), quote = F, row.names = F, col.names = F)
  write.table(x = down, file = sprintf('%s/%s_down.txt',out_path, md), quote = F, row.names = F, col.names = F)
  
  md_mets = md_mets_all
}

mets = md_mets
mets = mets[mets$q_value<=0.25,]
df = data.frame(metabolite = unique(mets$variable), module = '')
for ( i in 1:dim(df)[1] )
{
  temp = mets[mets$variable == df$metabolite[i],]
  df$module[i] = temp$module[temp$correlation == max(temp$correlation)]
}

