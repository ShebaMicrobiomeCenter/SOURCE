## filter metabolomics data using WGCNA modules
modules_mets_file = 'res_v2/SOURCE_Israel_TI_withStoolMetablolomicsV2LogClean_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default//selected_modules//SOURCE_Israel_TI_withStoolMetablolomicsV2LogClean_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default_module_cor_df.txt'
md_mets = read.table(file = modules_mets_file, header = T, sep = '\t')
md_mets$variable = make.names(md_mets$variable)
metadata_prms = c('Gender','Age_years','BMI','Disease','inflammed','CRP_numeric','Calprotectin_numeric','Active_Smoker_manual','Previous_or_Current_Tobacco_Use_manual','health_index_same_sample','health_index_stool')
md_mets = md_mets[! md_mets$variable %in% c(metadata_prms),]

out_path = 'WGCNA_mets_modules_lists/'
dir.create(out_path)

wanted_modules = c('yellow','green','red','pink', 'purple','tan','salmon','black','brown')
for (md in wanted_modules)
{
  up_mets = md_mets$variable[ md_mets$q_value <= 0.25 & md_mets$module == md & md_mets$correlation > 0]
  down_mets = md_mets$variable[ md_mets$q_value <= 0.25 & md_mets$module == md & md_mets$correlation < 0]
  write.table(x = up_mets, file = sprintf('%s/%s_up.txt', out_path, md), quote = F,sep = '\t', row.names = F, col.names = F)
  write.table(x = down_mets, file = sprintf('%s/%s_down.txt', out_path, md), quote = F,sep = '\t', row.names = F, col.names = F)
}