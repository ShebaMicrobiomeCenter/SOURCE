source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
source('/pita/users/tzipi/code/R_figs/maaslin2_funs.R')

supps_path = '../../supp_tables/halla/mets/'
# supps_path = '../../supp_tables/halla/'
pths = list.dirs(path = supps_path, full.names = TRUE, recursive = F)
pths = pths[grepl(pattern = 'vs',x = pths)]

for ( p in pths )
{
  f = sprintf('%s/halla_res_alpha0.25/all_associations.txt', p)
  df = read.table(f, header = T, sep = '\t')
  asv_num = gsub('.*_A','',df$Y_features)
  temp = tax_num_2_sv(asv_num)
  df$full_ASV = temp$Feature_ID
  df$full_taxonoomy = temp$Taxon
  
  out_name = gsub( 'all_associations', 'all_associations_ASV', f )
  write.table(x = df, file = out_name, quote = F, sep = '\t', row.names = F)
  
}

