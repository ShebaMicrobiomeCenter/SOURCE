# writing the neccesary text file for maaslin2 run
# fits taxa merged dataframe to 2 datafrmaes for write_maaslin2_taxa_map
write_maaslin2 = function(taxa, metadata_cols, maas_path, taxonomy_flag = F, taxonomy_file = '/pita/users/tzipi/projects/16S/all_merges/DB1-20_merge/res/16S_DB1-20_merged/taxonomy/metadata_v2_source.tsv')
{
  # seprating taxa to taxonomy and metadata
  # pos = grep('__',names(taxa))[1]-1
  pos = which(names(taxa) == 'Description')
  metadata = taxa[, c('SampleID',metadata_cols) ]
  tx = taxa[ ,c( 1,(pos+1):dim(taxa)[2] ) ]
  
  ASV = write_maaslin2_taxa_map(tx, metadata, maas_path, taxonomy_flag, taxonomy_file)
  return(ASV)
}

# write files neccesary for maaslin2 R package, using fitted taxonomy and metadata dataframes
# if wanted (taxonomy_flag) change the svs to taxonomy last 2 levels + fitted number
write_maaslin2_taxa_map = function(tx, metadata, maas_path, taxonomy_flag = F, taxonomy_file = '/pita/users/tzipi/projects/16S/all_merges/DB1-20_merge/res/16S_DB1-20_merged/taxonomy/metadata_v2_source.tsv')
{
  dir.create(maas_path)
  # writing to text files
  metadata_file = sprintf('%s/metadata.tsv',maas_path)
  write.table(metadata, metadata_file, sep = '\t',quote = F,row.names = F)
  tx_file = sprintf('%s/taxonomy.tsv',maas_path)
  if (!taxonomy_flag) 
  {
    write.table(tx, tx_file, sep = '\t',quote = F,row.names = F)
    svs = names(tx)[2:dim(tx)[2]]
    return(data.frame(ASV = svs))
  } else
  { # chamge svs to taxonomic annotation ( last 2 levels+number asociated to tx file full svs )
    svs = names(tx)[2:dim(tx)[2]]
    svs_df = sv_2_tax_v2(svs, taxonomy_file = taxonomy_file)
    # txn_names = sprintf('%s_A%s', clean_tax_names(svs_df$Taxon), svs_df$SV_number )
    
    txn = tx
    names(txn) = c('SampleID',svs_df$clean_name)
    write.table(txn, tx_file, sep = '\t',quote = F,row.names = F)
    
    # tx_names = sprintf('%s_ASV%d', svs, (2:dim(tx)[2])-1 )
    # names(tx) = c('SampleID',tx_names)
    tx_file_svs = sprintf('%s/taxonomy_svs.tsv',maas_path)
    write.table(tx, tx_file_svs, sep = '\t',quote = F,row.names = F)
    
    svs_tax = data.frame(Taxonomy = svs_df$clean_name, ASV = svs_df$Feature_ID)
    return(svs_tax)
  }
}

maaslin_add_asv_to_res = function(ASV, out_path)
{
  names(ASV)[1] = 'feature'
  # f = sprintf('%s/significant_results.tsv',out_path)
  # out_file =  sprintf('%s/significant_results_asv.tsv',out_path)
  f = sprintf('%s/significant_results.tsv',out_path)
  out_file =  sprintf('%s/significant_results_asv.tsv',out_path)
  df = read.table(file = f, sep = '\t',header = T)
  df2 = merge(df, ASV, by='feature',all.x=T)
  df2$ASV = gsub(pattern = '_.*',replacement = '',x = df2$ASV)
  df2 =df2[order(match(df2$feature,df$feature)),]
  # temp = sv_2_tax( df2$ASV )
  # df2$full_taxonomy = as.character(temp[,2])
  write.table(x = df2, file = out_file, quote = F, sep = '\t',row.names = F)
}

# fit svs to into annotation using taonomy table genrated by qiime2. 
sv_2_tax = function(svs, taxonomy_file = '/pita/users/tzipi/projects/16S/all_merges/DB1-14_merge/res/16S_DB1-14_merged/taxonomy/metadata.tsv')
{
  txn = read.table(file = taxonomy_file, header = T, sep = '\t')
  temp = data.frame(Feature.ID = as.character(svs))
  res = merge(x = temp, y = txn, by = 'Feature.ID', all.x = TRUE)
  return(res)
}

# fit svs to into annotation using taonomy table genrated by qiime2. 
sv_2_tax_v2 = function(svs, taxonomy_file = '/pita/users/tzipi/projects/16S/all_merges/DB1-19_merge/res/16S_DB1-19_merged/taxonomy/metadata_v2.tsv')
{
  txn = read.table(file = taxonomy_file, header = T, sep = '\t')
  temp = data.frame(Feature_ID = as.character(svs))
  res = merge(x = temp, y = txn, by = 'Feature_ID', all.x = TRUE, sort = F)
  res$Taxon = gsub(pattern = ';_',replacement = '; ', res$Taxon)
  res$clean_name = sprintf('%s_A%s', clean_tax_names(res$Taxon), res$SV_number )
  return(res)
}

# fit svs to into annotation using taonomy table genrated by qiime2. 
## SV num should be of format SV00001
tax_num_2_sv = function(tax_num_list, taxonomy_file = '/pita/users/tzipi/projects/16S/all_merges/DB1-19_merge/res/16S_DB1-19_merged/taxonomy/metadata_v2.tsv')
{
  txn = read.table(file = taxonomy_file, header = T, sep = '\t')
  temp = data.frame(SV_number = as.character(tax_num_list), order = 1:length(tax_num_list))
  res = merge(x = temp, y = txn, by = 'SV_number', all.x = TRUE)
  res = res[order(res$order),] ## order was added to keep original order
  res = res[,names(res)!='order']
  res$Taxon = gsub(pattern = ';_',replacement = '; ', res$Taxon)
  return(res)
}

# make full taxonomic annotation name into R readble name, with inly last 2 annotated levels.
clean_tax_names = function(txn)
{
  tx = make.names(txn)
  for (i in 1:7)
    tx = gsub(pattern = '[.][.][a-z]__$',replacement = '',x = tx)
  temp = stringr::str_split(string = tx,pattern = '[.][.]')
  for (i in 1:length(temp))
  {
    l = length(temp[[i]])
    if ( l == 1 )
    {
      tx[i] = temp[[i]]
    } else
    {
      tx[i] = sprintf('%s.%s',temp[[i]][l-1], temp[[i]][l] )
    }
  }
  return(tx)
}


maaslin2_plot_colored_scatters = function(taxa, x_var, col_var, head_path, maas_name, out_dir_name = 'my_figures', fdr = 0.25)
{
  in_file = sprintf('%s/%s/res/significant_results_asv.tsv',head_path, maas_name)
  df = read.table(in_file, header = T)
  
  df = df[df$metadata == x_var,]
  df = df[df$qval <= fdr,]
  
  for ( i in 1:length(df$feature) )
  {
    tx = df$feature[i]
    # sv = tax_num_2_sv(gsub('.*_ASV','SV',tx), taxonomy_file)
    # sv = sv[[2]]
    sv = df$ASV[i]
    
    p = ggplot(taxa, aes_string(x=x_var,y=sv, colour = col_var)) + geom_point(size=3) + theme_bw() + 
      # xlab('Diameter') + 
      ggtitle(sprintf('coefficient = %.3f, Q.value = %.4f', df$coef[i], df$qval[i])) + 
      ylab(gsub('_ASV','\nASV',tx))
    p = format_fig(fig = p, set_ttl_flag = F)
    
    out_path = sprintf('%s/%s/res/%s',head_path,maas_name, out_dir_name)
    dir.create(out_path)
    out_fig = sprintf('%s/%s_%s_%s_scatter_col_%s.tiff', out_path, i, x_var, tx, col_var )
    ggsave(out_fig,plot = p, device = 'tiff', width = 7,height = 5, compression  = 'lzw')
    
  }
}

maaslin2_set_data_for_heatmap_plot = function(taxa, head_path, maas_name, x_var, taxa_sig_var, group_var = '', 
        taxonomy_file = '/pita/users/tzipi/projects/16S/all_merges/DB1-20_merge/res/16S_DB1-20_merged/taxonomy/metadata_v2_source.tsv', 
        fdr = 0.25, taxonomy_flag = T)
{
  # browser()
  feature_prm = 'ASV'
  if (taxonomy_flag) {
    in_file = sprintf('%s/%s/res/significant_results_asv.tsv',head_path, maas_name)
  } else {
    in_file = sprintf('%s/%s/res/significant_results.tsv',head_path, maas_name)
    # feature_prm = 'feature'
  } 
  
  df = read.table(in_file, header = T)
  if ( !taxonomy_flag ) { df[[feature_prm]] = df$feature }
  df = df[df$metadata == taxa_sig_var,]
  df = df[df$qval <= fdr,]

  names(taxa) = make.names(names(taxa))
  if (!all(df[[feature_prm]] %in% names(taxa)))
  {
    print('Some siginificant features missing from the taxa dataframe!!')
    df = df[df[[feature_prm]] %in% names(taxa),]
  }
  
  library(reshape)
  vrs = unique( c('SampleID',x_var, group_var, taxa_sig_var) )
  vrs = vrs[vrs !='']
  
  taxa_m = reshape::melt(data = taxa, id.vars = vrs,measure.vars = make.names(df[[feature_prm]]), variable_name = feature_prm )
  
  if (taxonomy_flag)
  {
    taxa_m$taxa = sv_2_tax_v2(svs = taxa_m[[feature_prm]], taxonomy_file = taxonomy_file)$clean_name
  } else {
    taxa_m$taxa = taxa_m[[feature_prm]]
  }

  taxa_m$log10_RA = log10(taxa_m$value + 0.0001) 
  
  ## order by coefficient and the variable of interest
  taxa_m = merge(taxa_m, df[,c('ASV','coef')], by='ASV')
  temp = unique(taxa_m[, c('taxa','coef')])
  taxa_m$taxa = factor(taxa_m$taxa, levels = temp$taxa[order(temp$coef)])
  taxa_m[[taxa_sig_var]] = factor(taxa_m[[taxa_sig_var]], levels = unique(taxa_m[[taxa_sig_var]])[order(unique(taxa_m[[taxa_sig_var]] ))])
  taxa_m[[x_var]] = factor(taxa_m[[x_var]], levels = unique(taxa_m[[x_var]])[order(unique(taxa_m[[x_var]] ))])
  
  ## order taxa using clustering
  # library(reshape2)
  # data <- reshape2::dcast(taxa_m, SampleID ~ taxa, value.var="log10_RA")
  # row.names(data) = data$SampleID; data = data[, names(data) != 'SampleID']
  # ord <- hclust( dist( t(data), method = "euclidean"), method = "ward.D" )$order
  # taxa_m$taxa = factor(taxa_m$taxa, levels = names(data)[ord ])

  return(taxa_m)
}

maaslin2_healtmap_plot = function(taxa, head_path, maas_name, x_var, group_var = '', 
                                  taxonomy_file = '/pita/users/tzipi/projects/16S/all_merges/DB1-20_merge/res/16S_DB1-20_merged/taxonomy/metadata_v2_source.tsv', 
                                  out_dir_name = 'my_figures/', add_name = '', taxa_sig_var = '', fdr = 0.25, taxonomy_flag = T)
{
  if (taxa_sig_var == '')
    taxa_sig_var = x_var
  taxa_m = maaslin2_set_data_for_heatmap_plot(taxa, head_path, maas_name, x_var, taxa_sig_var, group_var = group_var, taxonomy_file = taxonomy_file, fdr,taxonomy_flag = taxonomy_flag)

  library(viridis)
  p_heatmap = ggplot(taxa_m, aes_string(x=x_var, y='taxa', fill = 'log10_RA')) + geom_tile() + 
    theme_classic() + 
    ylab('ASV') + 
    scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_fill_viridis() 
  
  if (group_var != '')
    p_heatmap = p_heatmap + facet_grid(as.formula(sprintf('.~%s',group_var)), scales = 'free', space = 'free') 
  out_path = sprintf('%s/%s/res/%s',head_path,maas_name, out_dir_name)
  dir.create(out_path)
  out_fig = sprintf('%s/%s_sigTaxa_%s_%s_heatmap%s_fdr%s.tiff', out_path, taxa_sig_var, x_var, group_var, add_name, fdr)
  ggsave(out_fig,plot = p_heatmap, device = 'tiff', width = 12,height = 2+length(unique(taxa_m$ASV))/8, compression  = 'lzw')
  
  return(p_heatmap)
}

