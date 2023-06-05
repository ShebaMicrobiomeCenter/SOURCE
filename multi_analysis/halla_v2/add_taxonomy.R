source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
source('/pita/users/tzipi/code/R_figs/maaslin2_funs.R')

in_path = 'res/full_data/'
out_path = in_path

fs = list.files(path = in_path)
fs_clean = gsub('_data.txt','',fs)

loc = gsub('.*_','',fs_clean)
tp = gsub('_.*','',fs_clean)

for ( i in 1:length(fs) )
{
  if (! grepl(pattern = '16S_taxa',x = fs[i]))
  if ( tp[i] == '16S' )
  {
    df = read.table(sprintf('%s/%s',in_path,fs[i]), header = T, sep = '\t')

    temp = sv_2_tax_v2(row.names(df), 
                       taxonomy_file = '/pita/users/tzipi/projects/16S/all_merges/DB1-20_merge/res/16S_DB1-20_merged/taxonomy/metadata_v2_source.tsv')

    row.names(df) = temp$clean_name
    out_name = gsub( '16S', '16S_taxa', fs[i] )
    write.table(x = df, file = sprintf('%s/%s', out_path, out_name), quote = F, sep = '\t', row.names = T)
  }
}