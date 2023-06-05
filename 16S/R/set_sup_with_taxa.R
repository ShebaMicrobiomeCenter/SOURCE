source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
source('/pita/users/tzipi/code/R_figs/maaslin2_funs.R')

in_path = '../../supp_tables/16S_bioms/'
out_path = in_path

fs = c('china_stool_amnonFlt_source_deblur_map_SVs.1.txt','israel_noRS_amnonFlt_source_deblur_map_SVs.1.txt')


for ( i in 1:length(fs) )
{
  df = read.table(sprintf('%s/%s',in_path,fs[i]), header = T, sep = '\t')
  pos = which(names(df) == 'Description')
  df = df[,c( 'sample_ID', names(df)[ (pos+1):dim(df)[2] ])]
  df = as.data.frame(t(df))
  temp = sv_2_tax_v2(row.names(df)[2:length(row.names(df))], 
                     taxonomy_file = '/pita/users/tzipi/projects/16S/all_merges/DB1-20_merge/res/16S_DB1-20_merged/taxonomy/metadata_v2_source.tsv')
  
  df$taxonomy = c('Taxonomy',temp$clean_name)
  l = dim(df)[2]
  df = df[,c(l,1:l-1)]
  df = as.data.frame(t(df))
  names(df)[1] = 'ASV'
  out_name = gsub( 'SVs.1', 'SVs_taxa', fs[i] )
  write.table(x = df, file = sprintf('%s/%s', out_path, out_name), quote = F, sep = '\t', row.names = F)
}