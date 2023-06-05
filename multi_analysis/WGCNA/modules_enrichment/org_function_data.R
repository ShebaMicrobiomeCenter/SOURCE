
get_good_enrichment_names_list = function(df, top_num = 30)
{
  df = df[df$Category %in% c('GO: Molecular Function','GO: Biological Process','GO: Cellular Component',
                                'Pathway','Coexpression Atlas','ToppCell Atlas'),]
  good_names = c()
  cats = unique(df$Category)
  for( i in 1:length(cats) )
  {
    temp = df[df$Category == cats[i],]
    temp = temp[base::order(temp$q.value.FDR.B.H, decreasing =F),]
    temp = temp[1:top_num,]
    good_names = c(good_names,temp$Name)
  }
  good_names = unique(good_names)
  return(good_names)
}

wanted_modules = c('yellow','green','red','pink','purple','tan','salmon','black','brown')

top_data_file_ptr = 'rnaSeq_main/modules_enrichment/data/SOURCE_Israel_TI_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default/enrichment_files/toppgene_%s_ens.txt'
top_data_file = sprintf(top_data_file_ptr,wanted_modules[1])
dff = read.table(top_data_file, sep = '\t',header = T, comment.char = '', quote="\"", na.strings = c('NA','na','',' '))
# good_names = c(good_names, get_good_enrichment_names_list(dff))
dff$module = wanted_modules[1]
for (i in 2:length(wanted_modules) )
{
  top_data_file = sprintf(top_data_file_ptr,wanted_modules[i])
  temp = read.table(top_data_file, sep = '\t',header = T, quote = '', comment.char = '')
  # good_names = c(good_names, get_good_enrichment_names_list(temp))
  temp$module = wanted_modules[i]
  dff = rbind(dff, temp)
}
# good_names = unique(good_names)
# dff = dff[dff$Name %in% good_names,]
# fdr_cut = 8
# fdr_cut = 6

dff = dff[dff$Category %in% c('GO: Molecular Function','GO: Biological Process','GO: Cellular Component',
                              'Pathway','Coexpression Atlas','ToppCell Atlas'),]


dff$FDR_BH = -1*log10(dff$q.value.FDR.B.H)
# dff$FDR_BH[dff$FDR_BH < fdr_cut ] = NA
# # dff$FDR_BH[dff$FDR_BH >20 ] = 20


interesting_modules = wanted_modules

temp  = dff[,c('Name','FDR_BH','module','Category','ID')]
# temp$Name = gsub(',','_',temp$Name)
temp = reshape2::dcast(temp, Name + Category + ID ~ module, value.var = 'FDR_BH', fun.aggregate=sum)
# row.names(temp) = temp$Name
# temp = temp[,-1]
temp[is.na(temp)]=0

# temp$Name = make.names(temp$Name)

out_file = 'rnaSeq_main/modules_enrichment/data/SOURCE_Israel_TI_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default/enrichment_files/toppgene_merged_ens.txt'
write.table(x = temp, file = out_file, quote = F, row.names = F, col.names = T, sep='\t')

