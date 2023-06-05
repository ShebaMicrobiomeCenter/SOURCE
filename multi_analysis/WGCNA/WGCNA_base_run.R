library(WGCNA)
library(ggplot2)
options(stringsAsFactors = FALSE);

source('/pita/users/tzipi/code/R_figs/WGCNA_funcs.R')

### set parameters
# cor_type = 'spearman'
cor_type = 'WGCNA_default'
network_type = 'signed hybrid'
# network_type = 'unsigned'
# network_type = 'signed'
plot_disease_related_modules_flag = T

added_data_type = 'withFFQv10'
# added_data_type = 'withStoolMetablolomicsV2LogClean'
# added_data_type = 'no'

# cohort = 'SOURCE'
# location = 'Stool'
# main_name = 'SOURCE_Stool'

cohort = 'SOURCE'
# country = 'Israel'
country = 'China'
location = 'TI'
oneSample = 'oneSample'
# main_name = 'SOURCE_Israel_TI'
main_name = sprintf('%s_%s_%s',cohort, country, location)

inf_flag = F
if ( inf_flag )
  main_name = sprintf('%s_inflamed', main_name)

age_match_flag = T
if ( age_match_flag )
  main_name = sprintf('%s_ageMatch', main_name)



path = 'rnaSeq_main/res_v2/'
if (main_name == 'SOURCE_Israel_TI' & oneSample == 'oneSample' & network_type == 'signed hybrid')
{
  cut_hight = 70000
  pwr = 12
  automatic_power_flag = F
  minModuleSize = 30
}
if (main_name == 'SOURCE_China_TI_ageMatch' & oneSample == 'oneSample' & network_type == 'signed hybrid')
{
  cut_hight = 50000
  pwr = 14
  automatic_power_flag = F
  minModuleSize = 30
}

net_name = sprintf('%s_%s_netType_%s_minMdSz_%s', oneSample, ifelse(automatic_power_flag, 'AutoPwr',sprintf('pwr_%s',pwr)), network_type, minModuleSize)
name = sprintf('%s_%s_%s_%s', main_name, added_data_type, net_name, cor_type)
name = gsub(' ','_',name)
name = gsub('_no_','_',name)
path = sprintf('%s/%s',path, name)
dir.create(path)

### input
# set needed tables
metadata_file = '../../rnaSeq/China_Israel/data/SOURCE_israel_china_rnaSeq_map_v2.txt'
metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'))
metadata_df$SampleID = make.names(metadata_df$SampleID)
metadata_df$Cohort = 'SOURCE'
metadata_df = metadata_df[metadata_df$Country == country & metadata_df$Location == location,]

if (inf_flag)
  metadata_df = metadata_df[metadata_df$inflammation_status == 'Inflamed' | metadata_df$Dx == 'healthy',]

if ( age_match_flag )
{
  metadata_df = metadata_df[metadata_df$Age_years >= 25 & metadata_df$Age_years <= 53,]
}

## adding Amnon's health index
health_index_file = '../../metadata/health_index_amnon.txt'
health_index = read.table(file = health_index_file, header = T,sep = '\t')
metadata_df$sample_ID = sprintf('%s_%s',metadata_df$pn_ID, metadata_df$location_inflamation)

metadata_df2 = merge(x = metadata_df, y=health_index, by = 'sample_ID', all.x = T, all.y = F)
metadata_df$health_index_same_sample = metadata_df2$health_index

metadata_df3 = merge(x = metadata_df, y=health_index[health_index$location == 'Stool',], by = 'pn_ID', all.x = T, all.y = F)
metadata_df$health_index_stool = metadata_df3$health_index

metadata_df$inflammation_status[metadata_df$Dx == 'healthy'] = 'Non-inflamed'

## adding manual smoking status
smoking_file = '../../metadata/alona_data/SOURCE_Israel_China_data_v11.txt'
smoking_df = read.table(file = smoking_file, header = T,sep = '\t')

metadata_df2 = merge(x = metadata_df, y=smoking_df, by = 'pn_ID', all.x = T, all.y = F)
metadata_df$Active_Smoker_manual = metadata_df2$X12.Active_Smoker_manual
metadata_df$Previous_or_Current_Tobacco_Use_manual = metadata_df2$X12.Previous_or_Current_Tobacco_Use_manual

data_file = sprintf('../../rnaSeq/China_Israel/res/SOURCE_%s_%s_geneFiltered_txi_res.txt', country, location)
ftr_df = read.table(file = data_file, header = T,sep = '\t', stringsAsFactors = F,row.names = 1)

## filter to one sapmle per patient - if there are 2, keep the inflammed sample.
if ( oneSample == 'oneSample' ) 
{
  good_pos = c()
  pns = unique(metadata_df$pn_ID)
  for (pn in pns)
  {
    pos = which(metadata_df$pn_ID == pn)
    if (length(pos) == 1)
      good_pos = c(good_pos, pos)
    if (length(pos) == 2)
    {
      p = which(metadata_df$pn_ID == pn & metadata_df$inflammation_status=='Inflamed')
      if(length(p) == 1)
      {
        good_pos = c(good_pos, p)
      } else
        print(sprintf('more than one inflammed sample for patient %s   problem!', pn))
    }
    if (length(pos) > 2)
      print(sprintf('more than 2 samples for patient %s   problem!', pn))
  }
  metadata_df = metadata_df[good_pos,]
}

# sync metadata and feature data frames
metadata_df = metadata_df[ metadata_df$SampleID %in% names(ftr_df), ]
metadata_df = metadata_df[order(metadata_df$SampleID),]
ftr_df = ftr_df[, metadata_df$SampleID]
ftr_df = as.data.frame(t(ftr_df))

### add and organize metadata
wanted_metadata =  c('Age_years','Gender','BMI',
                     'Active_Smoker_manual','Previous_or_Current_Tobacco_Use_manual',
                     'Disease','inflammed','CRP_numeric','Calprotectin_numeric', # 'CRP_log','Calprotectin_log' ,
                     'health_index_same_sample','health_index_stool')
metadata_prms = c('Gender','Age_years','BMI','Energy_Kcal_day',
                  'Active_Smoker_manual','Previous_or_Current_Tobacco_Use_manual',
                  'Disease','inflammed', 'CRP_numeric',# 'CRP_log',
                  'Calprotectin_numeric',# 'Calprotectin_log',
                  'health_index_stool','health_index_same_sample')
# wanted_metadata =  c('Age_years','Gender','BMI','Disease','inflammed','CRP_numeric','health_index_same_sample','health_index_stool' )
# metadata_prms = c('Gender','Age_years','BMI','Energy_Kcal_day','Disease','inflammed','CRP_numeric','health_index_stool','health_index_same_sample')
metadata_df$Disease = ifelse(metadata_df$Dx=="healthy", 0, 1)
metadata_df$inflammed = ifelse(metadata_df$inflammation_status=="Inflamed", 1, 0)
metadata_df$Gender = ifelse(metadata_df$Gender=="male",1, 0)
metadata_df$Age_years = as.numeric(metadata_df$Age_years)
metadata_df$Active_Smoker_manual = ifelse(metadata_df$Active_Smoker_manual=="Yes",1, 0)
metadata_df$Previous_or_Current_Tobacco_Use_manual = ifelse(metadata_df$Previous_or_Current_Tobacco_Use_manual=="Yes",1, 0)
metadata_df$CRP_log = log10(metadata_df$CRP_numeric)
metadata_df$Calprotectin_log = log10(metadata_df$Calprotectin_numeric)
# colnames(metadata_df)[colnames(metadata_df) == 'Age_years'] = 'Age (years)'

metadata_df_old = metadata_df
row.names(metadata_df) = metadata_df$SampleID
metadata_df = metadata_df[,wanted_metadata]

if (grepl('FFQ',added_data_type))
{
  if (added_data_type == 'withFFQKcalNorm')
  {
    ffq_file = '../../metadata/FFQ/FFQ_v9_kcal_norm.txt'
  } else if ( added_data_type == 'withFFQabsolute' )
  {
    ffq_file = '../../metadata/FFQ/FFQ_v9_absolute.txt'
  } else if ( grepl('withFFQv10',added_data_type) )
  {
    ffq_file = '../../metadata/FFQ/FFQ_v10.txt'
  } 
  ffq = read.table(ffq_file, header = T, sep = '\t', row.names = 1)
  if( grepl('Clean',added_data_type) )
  {
    for ( i in 1:dim(ffq)[2] )
    {
      # # change zeroes to 5th of the smallest value per metabolite
      # ffq[ffq[,i] == 0,i] = min(ffq[ffq[,i] != 0,i], na.rm = T)/5
      # change values over 4 SDs from the mean to the top/bottom value without it
      top_cut = mean(ffq[,i], na.rm = T) + 4*sd(ffq[,i], na.rm = T)
      ffq[ffq[,i] > top_cut & !is.na(ffq[,i]),i ] = max(ffq[ ffq[,i]<=top_cut,i ], na.rm = T)
      bottom_cut = mean(ffq[,i], na.rm = T) - 4*sd(ffq[,i], na.rm = T)
      ffq[ffq[,i] < bottom_cut & !is.na(ffq[,i]),i ] = min(ffq[ ffq[,i]>=bottom_cut,i ], na.rm = T)
    }
  }
  if( grepl('Log',added_data_type) )
    ffq = log(ffq + 0.00001)
    
  metadata_df2 = metadata_df
  metadata_df2$pn_ID = gsub('S','A',metadata_df_old$pn_ID)
  ffq$pn_ID = row.names(ffq)
  metadata_df2 = merge(metadata_df2, ffq, by = 'pn_ID', all.x = T)
  row.names(metadata_df2) = metadata_df2$pn_ID
  metadata_df2 = metadata_df2[gsub('S','A',metadata_df_old$pn_ID),]
  row.names(metadata_df2) = metadata_df_old$SampleID
  metadata_df = metadata_df2[,names(metadata_df2) != 'pn_ID']
} else if( grepl('Metablolomics',added_data_type) )
{
  # # set needed tables  
  if( grepl('Stool',added_data_type) )
  {
    if ( grepl('V2',added_data_type) )
    {
      met_metadata_file = '../../metabolomics/run2_oct22/data/SOURCE2022_Stool_metadata.txt'
      met_data_file = '../../metabolomics/run2_oct22/data/SOURCE2022_Stool_Data_MZM13_N2_filtered_v2.csv'
    } else
    {
      # met_metadata_file = '/pita/users/tzipi/projects/multiomics/nina_metabolomics/data/nina_10may22/Stool_metadata.txt'
      met_metadata_file = '/pita/users/tzipi/projects/multiomics/nina_metabolomics/data/nina_10may22/Stool_metadata_v2.txt'
      met_data_file = sprintf('/pita/users/tzipi/projects/multiomics/nina_metabolomics/data/nina_10may22/metabolomics/Stool_Data_MZM10_N2_filtered_calour.csv')
    }
  } else if ( grepl('Serum',added_data_type) )
  {
    if ( grepl('V2',added_data_type) )
    {
      met_metadata_file = '../../metabolomics/run2_oct22/data/SOURCE2022_Serum_metadata.txt'
      met_data_file = '../../metabolomics/run2_oct22/data/SOURCE2022_Serum_Data_MZM13_N2_filtered_v2.csv'
    } else
    {
      met_metadata_file = '/pita/users/tzipi/projects/multiomics/nina_metabolomics/data/nina_10may22/Serum_metadata_v2.txt'
      met_data_file = sprintf('/pita/users/tzipi/projects/multiomics/nina_metabolomics/data/nina_10may22/metabolomics/Serum_Data_MZM10_N2_filtered_calour.csv')
    }
  }
  met_metadata_df = read.table(file = met_metadata_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'))
  met = read.csv(file = met_data_file, header = T,row.names = 1)
  
  met_metadata_df$SampleID = make.names(met_metadata_df$SampleID)
  met_metadata_df =met_metadata_df[met_metadata_df$Cohort == 'SOURCE',]
  
  met_metadata_df = met_metadata_df[ met_metadata_df$SampleID %in% names(met), ]
  met_metadata_df = met_metadata_df[order(met_metadata_df$SampleID),]
  met = met[, make.names(sprintf('%s', met_metadata_df$SampleID))]
  names(met) = gsub('S','A',met_metadata_df$pn_ID)
  met = as.data.frame(t(met))
  
  if( grepl('Log',added_data_type) )
  {
    if( grepl('LogClean',added_data_type) )
    {
      for ( i in 1:dim(met)[2] )
      {
        # change zeroes to 5th of the smallest value per metabolite
        met[met[,i] == 0,i] = min(met[met[,i] != 0,i], na.rm = T)/5
        # change values over 4 SDs from the mean to the top/bottom value without it
        top_cut = mean(met[,i], na.rm = T) + 4*sd(met[,i], na.rm = T)
        met[met[,i] > top_cut,i ] = max(met[ met[,i]<=top_cut,i ], na.rm = T)
        bottom_cut = mean(met[,i], na.rm = T) - 4*sd(met[,i], na.rm = T)
        met[met[,i] < bottom_cut,i ] = min(met[ met[,i]>=bottom_cut,i ], na.rm = T)
      }
      # log transform
      met = log10(met)
    } else
      met = log10(met + 0.00001)
  }

  met$pn_ID = row.names(met)
  
  metadata_df$pn_ID = gsub('S','A',metadata_df_old$pn_ID)
  metadata_df2 = metadata_df
  metadata_df2 = merge(metadata_df2, met, by = 'pn_ID', all.x = T)
  row.names(metadata_df2) = metadata_df2$pn_ID
  metadata_df2 = metadata_df2[gsub('S','A',metadata_df_old$pn_ID),]
  row.names(metadata_df2) = metadata_df_old$SampleID
  metadata_df = metadata_df2[,names(metadata_df2) != 'pn_ID']
} else if( grepl('16S',added_data_type) )
{
  # # set needed tables  
  if( grepl('Stool',added_data_type) )
  {
    taxa_file = '../../16S/res_israel_stool_amnonFlt/biom/source_deblur_map_SVs.1.txt'
  } 
  taxa_full = read.table(file = taxa_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'))
  taxa = taxa_full[,(which(names(taxa_full)=='Description')+1):dim(taxa_full)[2]]
  row.names(taxa) = taxa_full$pn_ID
  # taxa_metadata = taxa_full[,1:which(names(taxa_full)=='Description')]
  
  source('/pita/users/tzipi/code/R_figs/maaslin2_funs.R')
  temp = sv_2_tax_v2(names(taxa), 
                     taxonomy_file = '/pita/users/tzipi/projects/16S/all_merges/DB1-20_merge/res/16S_DB1-20_merged/taxonomy/metadata_v2_source.tsv')
  names(taxa) = temp$clean_name 
  
  if( grepl('Log',added_data_type) )
  {
    if( grepl('LogClean',added_data_type) )
    {
      for ( i in 1:dim(taxa)[2] )
      {
        # change zeroes to 5th of the smallest value per metabolite
        taxa[taxa[,i] == 0,i] = min(taxa[taxa[,i] != 0,i], na.rm = T)/5
        # change values over 4 SDs from the mean to the top/bottom value without it
        top_cut = mean(taxa[,i], na.rm = T) + 4*sd(taxa[,i], na.rm = T)
        taxa[taxa[,i] > top_cut,i ] = max(taxa[ taxa[,i]<=top_cut,i ], na.rm = T)
        bottom_cut = mean(taxa[,i], na.rm = T) - 4*sd(taxa[,i], na.rm = T)
        taxa[taxa[,i] < bottom_cut,i ] = min(taxa[ taxa[,i]>=bottom_cut,i ], na.rm = T)
      }
      # log transform
      taxa = log10(taxa)
    } else
      taxa = log10(taxa + 0.00001)
  }
  
  taxa$pn_ID = row.names(taxa)
  
  metadata_df$pn_ID = gsub('S','A',metadata_df_old$pn_ID)
  metadata_df2 = metadata_df
  metadata_df2 = merge(metadata_df2, taxa, by = 'pn_ID', all.x = T)
  row.names(metadata_df2) = metadata_df2$pn_ID
  metadata_df2 = metadata_df2[gsub('S','A',metadata_df_old$pn_ID),]
  row.names(metadata_df2) = metadata_df_old$SampleID
  metadata_df = metadata_df2[,names(metadata_df2) != 'pn_ID']
} else if( grepl('MTG',added_data_type) )
{
  temp = gsub('.*MTG','',added_data_type)
  temp_country = gsub('China','chinese',country) 
  temp_country = gsub('Israel','israeli',temp_country)
  ftr_file = sprintf('../../metagenomics/data/source_%s_mtg_%s.tsv',temp_country,temp)
  ftr = read.table(file = ftr_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'), comment.char = '', quote = '', row.names = 1)
  ftr = as.data.frame(t(ftr))
  
  ftr$pn_ID = row.names(ftr)
  
  metadata_df$pn_ID = gsub('S','A',metadata_df_old$pn_ID)
  metadata_df2 = metadata_df
  metadata_df2 = merge(metadata_df2, ftr, by = 'pn_ID', all.x = T)
  row.names(metadata_df2) = metadata_df2$pn_ID
  metadata_df2 = metadata_df2[gsub('S','A',metadata_df_old$pn_ID),]
  row.names(metadata_df2) = metadata_df_old$SampleID
  metadata_df = metadata_df2[,names(metadata_df2) != 'pn_ID']
}



## remove metadata columns with no data or no variance
# for (i in dim(metadata_df)[2]:1)
#   if ( sum(!is.na(metadata_df[,i])) <= 1 | var(metadata_df[,i], na.rm = T) ==0 ) 
#   { print(sprintf('remove %s from metadata',names(metadata_df)[i]))
#     metadata_df = metadata_df[,-i] }
metadata_df = clean_metadata_eampty_noVar(metadata_df)
                                       
### run QC and filter bad samples and features
ftr_df = WGCNA_initial_QC(ftr_df)
keepSamples = WGCNA_sample_tree(ftr_df, cut_hight, path, name)
ftr_df = ftr_df[keepSamples, ]
bad_samples = row.names(metadata_df)[!keepSamples]
metadata_df = metadata_df[row.names(ftr_df), ]

## plot new clean tree with metadata
WGCNA_sample_tree_with_metadata(ftr_df, metadata_df, path, name)

#### network construction
# Choose a set of soft-thresholding powers
sft = WGCNA_power_options_plot(ftr_df)
res = WGCNA_net(ftr_df, sft, automatic_power_flag, pwr, network_type, minModuleSize)
net = res[[1]]
# name = sprintf('%s_%s',name,res[[2]])
WGCNA_plot_net_dendro(net, path, name)
# table(net$colors)
moduleColors = WGCNA::labels2colors(net$colors)

## plot module correlations
# WGCNA_plot_genes_correlaitons_by_module(ftr_df, moduleColors, path)

## ME calculation
# Calculate MEs with color labels
res = WGCNA_calculate_module_eigengenes_and_corrs(ftr_df,moduleColors, metadata_df, cor_type, name)
MEs = res[[1]]
moduleTraitCor = res[[2]]
moduleTraitPvalue = res[[3]]


## plot ME/metdata corr heatmap
# original version
WGCNA_ME_met_heatmap_original(moduleTraitCor, moduleTraitPvalue, name, path)
if (plot_disease_related_modules_flag)
  WGCNA_ME_met_heatmap_original_prm_related(moduleTraitCor, moduleTraitPvalue, name, prm = 'Disease', cutoff = 0.05, path)

# my version
heatmap = WGCNA_ME_met_heatmap(moduleTraitCor, moduleTraitPvalue, name, metadata_prms)
ggsave(sprintf('%s/%s_module_cor_heatmap.pdf',path, name),plot = heatmap, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+dim(moduleTraitPvalue)[2]/10, limitsize = F)
if (plot_disease_related_modules_flag)
{
  heatmap2 = WGCNA_ME_met_heatmap_prm_related(moduleTraitCor, moduleTraitPvalue, name, prm = 'Disease', cutoff = 0.05, metadata_prms)
  ggsave(sprintf('%s/%s_module_cor_heatmap_%sF%s.pdf',path, name,'Disease',0.05), plot = heatmap2, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+dim(moduleTraitPvalue)[2]/10, limitsize = F)
  # heatmap3 = WGCNA_ME_met_heatmap_prm_related_pflt(moduleTraitCor, moduleTraitPvalue, name, prm = 'Disease', cutoff = 0.05, metadata_prms, p_cutoff = 0.0001)
  # y_lables = ggplot_build(heatmap3)$layout$panel_params[[1]]$y$get_labels()
  # ggsave(sprintf('%s/%s_module_cor_heatmap_%sF%s_pF%s.pdf',path, name,'Disease',0.05, 0.0001), plot = heatmap3, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+length(y_lables)/10, limitsize = F)
  res = WGCNA_ME_met_heatmap_full_res(moduleTraitCor,
                                      moduleTraitPvalue, name, metadata_prms, 
                                      min_p = 0.25, 
                                      fdr_flag_non_metadata = T, show_p_value_fdr_flag = T )
  heatmap4 = res[[1]]
  y_lables = ggplot_build(heatmap4)$layout$panel_params[[1]]$y$get_labels()
  ggsave(sprintf('%s/%s_module_cor_heatmap_q%s.pdf',path, name, 0.25),plot = heatmap4, device = 'pdf', width = 5+dim(moduleTraitPvalue)[1]/2,height = 5+length(y_lables)/10, limitsize = F)
  
}

## write output table
df = WGCNA_set_results_table(ftr_df, metadata_df, MEs)
write.table(x = df, file = sprintf('%s/%s_gene_module_table.txt',path, name),quote = F,row.names = F,sep = '\t')

## save relevant data frames to check things later
save(ftr_df, metadata_df, moduleColors, cor_type, name, MEs, moduleTraitCor, moduleTraitPvalue, heatmap, heatmap2, df, file = sprintf('%s/data.RData',path))
