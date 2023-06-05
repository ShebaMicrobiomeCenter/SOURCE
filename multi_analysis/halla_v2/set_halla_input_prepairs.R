source('/pita/users/tzipi/code/R_figs/figs_funcs.R')

taxonomy_file = '/pita/users/tzipi/projects/16S/all_merges/DB1-20_merge/res/16S_DB1-20_merged/taxonomy/metadata_v2_source.tsv'

# countries = c('Israel','China')
# country = 'Israel'
# 
# locations = c('TI','R','Stool','Serum')
# location = 'Stool'
# 
# types = c('16S','metabolomics','FFQ') #'rnaSeq'
# type = 'metabolomics'

get_exp_map = function(type, location = '', country ='', sv_to_taxa_fag = T)
{
  na_str = c('na','NA')
  if (type == 'metabolomics' & country == 'Israel' & location == 'Serum' ) ## using old data, need to fit to Nina's new analysis
  {
    ## metabolomics
    met_file = '../../metabolomics/run2_oct22/data/SOURCE2022_Serum_Data_MZM13_N2_filtered_v2.csv' 
    met = read.csv(file = met_file, header = T,row.names = 1)
    row.names(met) = make.names(row.names(met))
    
    map_met_file = '../../metabolomics/run2_oct22/data/SOURCE2022_Serum_metadata.txt'
    map_met = read.table(map_met_file, header = T, sep = '\t', row.names = 1, na.strings = na_str)
    map_met = map_met[map_met$Cohort == 'SOURCE',]
    row.names(map_met) = make.names(row.names(map_met))
    
    map_met = map_met[row.names(map_met) %in% names(met),]
    met_f = met[, row.names(map_met)]
    
    map_met$pn_ID = gsub('S','A',map_met$pn_ID)
    names(met_f) = map_met$pn_ID
    
    exp_df = met_f
    map_df = map_met
    name = 'Israel_metabolomics_Serum'
  }
  else if (type == 'metabolomics' & country == 'Israel' & location == 'Stool' ) 
  {
    ## metabolomics
    met_file = '../../metabolomics/run2_oct22/data/SOURCE2022_Stool_Data_MZM13_N2_filtered_v2.csv'
    met = read.csv(file = met_file, header = T,row.names = 1)
    row.names(met) = make.names(row.names(met))
    
    map_met_file = '../../metabolomics/run2_oct22/data/SOURCE2022_Stool_metadata.txt'
    map_met = read.table(map_met_file, header = T, sep = '\t', row.names = 1, na.strings = na_str)
    map_met = map_met[map_met$Cohort == 'SOURCE',]
    row.names(map_met) = make.names(row.names(map_met))
    
    map_met = map_met[row.names(map_met) %in% names(met),]
    met_f = met[, row.names(map_met)]
    
    map_met$pn_ID = gsub('S','A',map_met$pn_ID)
    names(met_f) = map_met$pn_ID
    
    exp_df = met_f
    map_df = map_met
    name = 'Israel_metabolomics_Stool'
  } 
  else if (type == 'metabolomics' & country == 'China' & location == 'Stool' ) 
  {
    ## metabolomics
    met_file = '../../metabolomics/China/data/SOURCE_china_mets_table_norm_v2.tsv'
    met = read.table(file = met_file, header = T,row.names = 1, sep='\t', quote = '', comment.char = '')
    row.names(met) = make.names(row.names(met))
    
    map_met_file = '../../metabolomics/China/data/SOURCE_china_mets_map.tsv'
    map_met = read.table(map_met_file, header = T, sep = '\t', row.names = 1, na.strings = na_str)
    row.names(map_met) = make.names(row.names(map_met))
    
    map_met = map_met[row.names(map_met) %in% names(met),]
    met_f = met[, row.names(map_met)]
    
    # map_met$pn_ID = gsub('S','A',map_met$pn_ID)
    # names(met_f) = map_met$pn_ID
    
    exp_df = met_f
    map_df = map_met
    name = 'China_metabolomics_Stool'
  } 
  else if ( type == 'rnaSeq' )
  {
    ## rnaSeq
    map_tpm_file = '../../rnaSeq/China_Israel/data/SOURCE_israel_china_rnaSeq_map_v2.txt'
    map_tpm = read.table(map_tpm_file, header = T, sep = '\t', row.names = 1, na.strings = na_str)
    
    ## filter to wanted country, location 
    map_tpm = map_tpm[grepl(pattern = country, x = map_tpm$Country),]
    map_tpm = map_tpm[map_tpm$Location == location,]
    
    tpm_file = sprintf('../../rnaSeq/China_Israel/res/SOURCE_%s_%s_geneFiltered_txi_res.txt',country,location)
    tpm = read.table(tpm_file, header = T, sep = '\t', row.names = 1, na.strings = na_str, comment.char = '', quote="\"")
    
    tpm = tpm[,make.names(row.names(map_tpm))]
    row.names(map_tpm) = sprintf('%s_%s',map_tpm$pn_ID, map_tpm$location_inflamation)
    names(tpm) = row.names(map_tpm)
    
    # temp = names(tpm)
    # names(tpm) = sprintf('%s_%s',substr(temp, 1, 4)  ,substr(temp, 5, nchar(temp))  )
    
    # DE_list_file = 'data/rnaSeq_v_lncRNA_TI/DE_gene_list_yael_8nov21.txt'
    # DE_list = read.table(DE_list_file, header = F, sep = '\t')
    # DE_list$V1 = DE_list[ ! DE_list %in% c('GUC','CSF2'),]
    
    map_tpm = map_tpm[sort(row.names(map_tpm)),]
    tpm = tpm[,row.names(map_tpm)]
    tpm_f = tpm
    
    # source('set_rnaSeq_module_PC1.R') ## adding wgcna PC values
    # # dx_modules = c('tan','yellow','brown','green','pink')
    # # dx_modules = c('blue','green','brown','salmon','tan')
    # modules_file = sprintf('../../rnaSeq/China_Israel/res/WGCNA_v2/SOURCE_%s_%s_all_oneSample_gene_module_table.txt', country, location)
    # mds = get_modules_PC1(tpm=tpm, modules_file = modules_file)
    # 
    # # tpm_f = tpm[DE_list$V1,]
    # tpm_f = rbind(mds, tpm)
    # # names(tpm_f) = map_tpm$pn_ID
    
    exp_df = tpm_f
    map_df = map_tpm
    name = sprintf('%s_%s_%s', country, type, location = '')
  } 
  else if ( type == '16S' )
  {
    if ( country == 'Israel' )
      taxa_file = '../../16S/res_israel_all_amnonFlt/biom/source_deblur_map_SVs.1.txt'
    if ( country == 'China' )
      taxa_file = '../../16S/res_china_stool_amnonFlt/biom/source_deblur_map_SVs.1.txt'
    
    na_str = c('no_data','_','NA','unknown', 'other','na','No_followup','ND')
    taxa = read.table(taxa_file, header = T, sep = '\t', na.strings = na_str, comment.char = '', quote="\"")
    
    ## filter to wanted country, location 
    # taxa = taxa[grepl(pattern = country, x = taxa$Cohort),]
    taxa = taxa[taxa$location == location,]
    
    map_taxa = taxa[,1:which(names(taxa)=='Description')]
    taxa_f = taxa[ ,(which(names(taxa)=='Description')+1):dim(taxa)[2]]
    # taxa_f = taxa[,DE_list$V1]
    
    taxa_f = as.data.frame(t(taxa_f))
    if ( location == 'Stool' )
    {
      names(taxa_f) = sprintf('%s',map_taxa$pn_ID)
    } else
      names(taxa_f) = sprintf('%s_%s',map_taxa$pn_ID, map_taxa$Location_inflamation)
    
    if (sv_to_taxa_fag)
    {
      source('/pita/users/tzipi/code/R_figs/maaslin2_funs.R')
      svs = row.names(taxa_f)
      svs_df = sv_2_tax_v2(svs, taxonomy_file = taxonomy_file)
      row.names(taxa_f) = svs_df$clean_name
    }
  
    exp_df = taxa_f
    map_df = map_taxa
    name = sprintf('%s_%s_%s', country, type, location)
  }
  else if ( grepl('metagenomics',type) )
  {
    if ( country == 'Israel' )
      taxa_file = sprintf('../../metagenomics/data/source_israeli_mtg_%s.tsv',gsub('.*_','',type))
    if ( country == 'China' )
      taxa_file = sprintf('../../metagenomics/data/source_chinese_mtg_%s.tsv',gsub('.*_','',type))
    
    na_str = c('no_data','_','NA','unknown', 'other','na','No_followup','ND')
    df = read.table(taxa_file, header = T, sep = '\t', na.strings = na_str, comment.char = '', quote="\"", row.names = 1)
    
    metadata_file = '../../metagenomics/data/SOURCE_Israel_China_data_v12.tsv'
    map = read.table(metadata_file,header = T,sep = '\t')
    
    row.names(map) = map$mtg_ID
    map = map[names(df),]
    
    names(df) = map$pn_ID
    
    # df_f = as.data.frame(t(df))
    df_f = df
    
    exp_df = df_f
    map_df = map
    name = sprintf('%s_%s', country, type)
  }
  else if (type == 'FFQ')
  {
    ffq_file = '../../metadata/FFQ/FFQ_v11.txt'
    ffq = read.table(ffq_file, header = T, sep = '\t', row.names = 1)
    # remove the questioner data, as it not real continuous.
    ffq = ffq[, !grepl('^X15',names(ffq))]
    
    metadata_file = '../../metadata/alona_data/SOURCE_Israel_China_data_v11.txt'
    md = read.table(metadata_file, header = T, sep='\t')
    
    md_country = ifelse( grepl('Israel',md$Patient_group), 'Israel','China' )
    md = md[md_country == country, ]
    md = md[,1:5]
    
    ffq = ffq[row.names(ffq) %in% md$pn_ID, ]
    
    source('/pita/users/tzipi/code/R_figs/WGCNA_funcs.R')
    ffq = clean_metadata_eampty_noVar(ffq)

    ffq = as.data.frame(t(ffq))
    md = md[md$pn_ID %in% names(ffq),]
    
    exp_df = ffq
    map_df = md
    name = sprintf('%s_%s', country, type)
  }
  if (!exists('exp_df'))
  {
    print(sprintf('No data for %s %s %s', type, location, country))
    return()
  }
  return(list(exp_df,map_df,name))
}

# res = get_exp_map(type, location, country)
# exp = res[[1]]
# map = res[[2]]
# name = res[[3]]

# out_path = 'res/prepairs_data/'
# write.table(x = exp, file = sprintf('%s/%s_data.txt', out_path, name), quote = F, sep = '\t', row.names = T)
