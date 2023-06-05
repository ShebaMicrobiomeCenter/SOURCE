
pth = 'rnaSeq_main/modules_enrichment/data/SOURCE_Israel_TI_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default/'
md_file = sprintf('%s/SOURCE_Israel_TI_withFFQv10_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default_gene_module_table.txt',pth)
md = read.table(md_file, header = T)

wanted_modules = c('yellow','green','red','pink','purple','tan','salmon','black','brown')

for ( m in wanted_modules )
{
  df = data.frame(Gene = md$Gene[md$module_color == m])
  df$gene_name = gsub('_.*','',df$Gene)
  df$ens_id = gsub('.*_','',df$Gene)
  write.table(x = df, file = sprintf('%s/genes_files/%s',pth,m), quote = F, row.names = F, col.names = T, sep='\t')
}

