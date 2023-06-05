
ftr_file = '../data/SOURCE_china_mets_table_t.tsv'

na_str = c('no_data','_','NA','unknown', 'other','na','No_followup')
ftr = read.table(file = ftr_file, header = T,row.names = 1, na.strings = na_str, sep='\t')

ftr_old = ftr
for ( i in 1:dim(ftr)[2] )
{
  ftr[,i] = ftr[,i]/sum(ftr[,i])*100
}


write.table(x = ftr, file = '../data/SOURCE_china_mets_table_t_norm.tsv', quote = F, sep = '\t', row.names = T)
