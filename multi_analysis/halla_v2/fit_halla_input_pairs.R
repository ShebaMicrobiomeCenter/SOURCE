source('/pita/users/tzipi/code/R_figs/figs_funcs.R')

## to handle one inflamed and one non-inflamed sample per patient
## the sample names must be of the shape "A000_TI1"
get_one_samples_per_patient = function(good_samples)
{
  good_pos = c()
  pn_ID = gsub(pattern = '_.*',replacement = '',good_samples)
  inflammation_status = ifelse(test = grepl(pattern = 'R3|TI1',x = good_samples),'Inflamed','non-inflamed')
  pns = unique( pn_ID )
  for (pn in pns)
  {
    pos = which(pn_ID == pn)
    if (length(pos) == 1)
      good_pos = c(good_pos, pos)
    if (length(pos) == 2)
    {
      p = which(pn_ID == pn & inflammation_status == 'Inflamed')
      if(length(p) == 1)
      {
        good_pos = c(good_pos, p)
      } else
        print(sprintf('more than one inflammed sample for patient %s   problem!', pn))
    }
    if (length(pos) > 2)
      print(sprintf('more than 2 samples for patient %s   problem!', pn))
  }
  return( good_samples[good_pos] )
}

countries = c('Israel','China')
locations = c('TI','R','Stool','Serum')
types = c('16S','metabolomics','FFQ','metagenomics_species','metagenomics_pathabundance','metagenomics_ecs') #'rnaSeq'

source('set_halla_input_prepairs.R')
       
country = 'Israel'
# country = 'China'

# res1 = get_exp_map('metabolomics', 'Stool', country)
# res1 = get_exp_map('metagenomics_pathabundance', 'Stool', country)
res1 = get_exp_map('metabolomics', 'Stool', country)
# res1 = get_exp_map('metabolomics', 'Serum', country)
# res1 = get_exp_map('FFQ', 'Stool', country, sv_to_taxa_fag = F)
# res1 = get_exp_map('16S', 'Stool', country)
d1 = res1[[1]]
map1 = res1[[2]]
name1 = res1[[3]]

# res2 = get_exp_map('16S', 'Stool', country, sv_to_taxa_fag = F)
# res2 = get_exp_map('FFQ', 'Stool', country, sv_to_taxa_fag = F)
# res2 = get_exp_map('metagenomics_ecsNamed', 'Stool', country)
res2 = get_exp_map('metagenomics_koNamed', 'Stool', country)
# res2 = get_exp_map('metagenomics_koL3', 'Stool', country)
d2 = res2[[1]]
map2 = res2[[2]]
# name2 = sprintf('%s_sv',res2[[3]])
name2 = res2[[3]]

# # res2 = get_exp_map('16S', 'Stool', country)
# res2 = get_exp_map('FFQ', 'Stool', country)
# d2 = res2[[1]]
# map2 = res2[[2]]
# name2 = res2[[3]]


## filter metabolomics data using WGCNA modules
modules_mets_file = '../WGCNA/rnaSeq_main/res_v2/SOURCE_Israel_TI_withStoolMetablolomicsV2LogClean_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default//selected_modules//SOURCE_Israel_TI_withStoolMetablolomicsV2LogClean_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default_module_cor_df.txt'
md_mets = read.table(file = modules_mets_file, header = T, sep = '\t')
md_mets$variable = make.names(md_mets$variable)
metadata_prms = c('Gender','Age_years','BMI','Disease','inflammed','CRP_numeric','Calprotectin_numeric','Active_Smoker_manual','Previous_or_Current_Tobacco_Use_manual','health_index_same_sample','health_index_stool')
md_mets = md_mets[! md_mets$variable %in% c(metadata_prms),]

md_mets_all = md_mets
d1_all = d1
name1_all = name1


# md_mets = md_mets[md_mets$q_value <= 0.25,]
# d1 = d1[unique(md_mets$variable),]
# name1 = gsub(pattern = 'metabolomics','metabolomicsAllModules',name1)

wanted_modules = c('yellow','green','red','pink', 'purple','tan','salmon','black','brown')
# for (md in wanted_modules)
{
  d1 = d1_all
  name1 = name1_all
  md_mets = md_mets_all

  # md_mets = md_mets[md_mets$q_value <= 0.25 & md_mets$module == md,]
  # d1 = d1[unique(md_mets$variable),]
  # name1 = gsub(pattern = 'metabolomics',sprintf('metabolomicsM%s',md),name1)

  
  # in_path = 'res/prepairs_data/'
  out_path = 'res/full_pairs/'
  dir.create(out_path)
  
  
  if ( all(grepl(pattern = '_',x = names(d1))) && all(grepl(pattern = '_',x = names(d2))) ) ## both data sets have more than one sample per patient
  {
    if ( loc[i] == loc[j] ) # same location, can directly use the same samples
    {
      good_samples = intersect(names(d1),names(d2))
      good_samples = get_one_samples_per_patient(good_samples)
    } else ## different locations, need to filter to one sample per patient and fit by pn_ID
    {
      d1 = d1[, get_one_samples_per_patient(names(d1)) ]
      d2 = d2[, get_one_samples_per_patient(names(d2)) ]
      names(d1) = gsub('_.*','',names(d1))
      names(d2) = gsub('_.*','',names(d2))
      good_samples = intersect(names(d1),names(d2))
    }
    
  } else if ( all(grepl(pattern = '_',x = names(d1))) ) ## only d1 has more than one sample per patient
  {
    temp = get_one_samples_per_patient(good_samples = names(d1))
    d1 = d1[,temp]
    names(d1) = gsub(pattern = '_.*','',names(d1))
    good_samples = intersect(names(d1),names(d2))
  } else if ( all(grepl(pattern = '_',x = names(d2))) ) ## only d2 has more than one sample per patient
  {
    temp = get_one_samples_per_patient(names(d2))
    d2 = d2[,temp]
    names(d2) = gsub(pattern = '_.*','',names(d2))
    good_samples = intersect(names(d1),names(d2))
  } else
  {
    good_samples = intersect(names(d1),names(d2))
  }
  ## filtering to one samples per patient
  d1 = d1[,good_samples]
  d2 = d2[,good_samples]
  
  if (dim(d1)[2] > 0 )
  {
    name = sprintf('%s_vs_%s',name1,name2)
    dir.create(sprintf('%s/%s/',out_path, name))
    
    write.table(x = d1, file = sprintf('%s/%s/x.txt', out_path,name), quote = F, sep = '\t')
    write.table(x = d2, file = sprintf('%s/%s/y.txt', out_path,name), quote = F, sep = '\t')
  }
  
}
