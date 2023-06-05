source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
source('/pita/users/tzipi/code/R_figs/maaslin2_funs.R')

supps_path = '../../supp_tables/halla/mets/'
# supps_path = '../../supp_tables/halla/'
pths = list.dirs(path = supps_path, full.names = TRUE, recursive = F)
pths = pths[grepl(pattern = 'vs',x = pths)]

types = c('16S_Stool','FFQ',
          'metagenomics_ecsNamed','metagenomics_pathabundance','metagenomics_species')
new_types_names = c('16S_Stool_taxa','FFQ',
                    'MGX_enzymes','MGX_pathways','MGX_taxa')

dfs = list()
for (j in 1:length(types))
{
  pths_f = pths[grepl(types[j], pths)]
  df_m = data.frame()
  for ( i in 1:length(pths_f) )
  {
    if (grepl( '16S',pths_f[i] ) )
    {
      f = sprintf('%s/halla_res_alpha0.25/all_associations_ASV.txt', pths_f[i])
    } else
      f = sprintf('%s/halla_res_alpha0.25/all_associations.txt', pths_f[i])
    df = read.table(f, header = T, sep = '\t', comment.char = '', quote = '')
    temp = gsub('.*metabolomicsM','',f); temp = gsub('_.*','',temp)
    df$metabolites_module = temp
    df$type = new_types_names[j]
    df_m = rbind(df_m, df)
  }
  dfs[[j]] = df_m
  write.table(x = df_m, file = sprintf('%s/halla_%s.txt', supps_path, new_types_names[j]), 
              row.names = F, quote = F, sep = '\t')
  write.table(x = df_m[df_m$q.values<=0.25,], file = sprintf('%s/halla_%s_q25.txt', supps_path, new_types_names[j]), 
              row.names = F, quote = F, sep = '\t')
}



# library(xlsx)
# eg <- list("one" = data.frame(one = rep(1, 100)),
#            "two" = data.frame(two = rep(2, 200)))
# wb <- createWorkbook()
# for (i in seq_along(eg)) {
#   sheet <- createSheet(wb, names(eg)[i])
#   addDataFrame(eg[i], sheet)
# }
# saveWorkbook(wb, "eg.xlsx")