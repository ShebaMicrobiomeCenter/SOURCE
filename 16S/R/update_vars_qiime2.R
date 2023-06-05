
library(ggplot2)

source('/pita/users/tzipi/code/R_figs/figs_funcs.R')


name = 'source'
# path = '../res_china_stool_amnonFlt/'
# path = '../res_israel_stool_amnonFlt/'
# path = '../res_israel_all_amnonFlt/'
# path = '../res_israel_noRS_amnonFlt/'
# path = '../res_israel_china_all_amnonFlt/'
path = '../res_china_stool_rural_ageMatched_amnonFlt/'


taxa_file = sprintf('%s/taxa-bar-plots/level-%%d.csv',path)
a_file = sprintf('%s/core-metrics-results/faith_pd_vector/alpha-diversity.tsv',path)
b_file = sprintf('%s/core-metrics-results/unweighted_unifrac_distance_matrix/distance-matrix.tsv',path)
pcoa_file = sprintf('%s/core-metrics-results/unweighted_unifrac_pcoa_results/ordination.txt',path)
out_file = sprintf('%s/biom/%s_deblur_map_L7.1.txt',path, name)


read_taxa_level = 7

read_file = sprintf(taxa_file, 2)


#read files
na_str = c('no_data','_','NA','unknown', 'other','na')
temp = read.table(read_file,sep=",", header=TRUE, na.strings = na_str, comment.char = '', quote="\"")
names(temp)[1] = 'SampleID'

# taxa_pos = grep(pattern = 'k__',names(temp))
# metadta_pos = grep(pattern = '^((?!__).)*$',names(temp),perl = T)

taxa_pos = grep(pattern = 'k__',names(temp))
metadata_pos = which ( !names(temp) %in% names(temp)[taxa_pos]  )


taxa = temp[,c(1,taxa_pos)]

abs_table_to_relative = function(taxa)
{
  good_pos = which( names(taxa) != 'k__Bacteria.__' & names(taxa) != 'k__Bacteria.__.__')
  taxa = taxa[,good_pos]
  
  for ( i in 1:dim(taxa)[1] )
    taxa[i,2:dim(taxa)[2]] = taxa[i,2:dim(taxa)[2]]/sum(taxa[i,2:dim(taxa)[2]])
  
  return(taxa)
}

# make relative abundance table
taxa = abs_table_to_relative(taxa)

metadata = temp[,metadata_pos]



for (i in 3:read_taxa_level)
{
  read_file = sprintf(taxa_file, i)
  temp = read.table(read_file,sep=",", header=TRUE, na.strings = na_str, comment.char = '', quote="\"")
  taxa_pos = grep(pattern = 'k__',names(temp))
  
  temp = temp[,c(1,taxa_pos)]
  names(temp)[1] = 'SampleID'
  
  # make relative abundance table
  temp = abs_table_to_relative(temp)
  
  taxa = merge(taxa, temp, by='SampleID')
}
taxa = merge(metadata, taxa, by='SampleID')


alphas = c('faith')

# add alpha diversity
for ( alph in alphas )
{
  read_file = a_file
  temp = read.table(read_file,sep="\t", header=TRUE, na.strings = na_str)
  names(temp)[1] = 'SampleID'

  taxa = merge(temp, taxa, by='SampleID')
}

# add dysbiosis index
mdi = calc_MD_index(taxa)

temp = data.frame(SampleID = taxa$SampleID, dysbiosis_index = mdi)
taxa = merge(temp, taxa, by='SampleID')


# resp = vector(mode = 'character',length = dim(taxa)[1])
# resp[taxa$response == 0] = 'noR'
# resp[taxa$response == 1] = 'R'
# resp[taxa$response == 9] = 'not_determined'
# resp[taxa$response == 'Control'] = 'Control'
# resp[is.na(taxa$response)] = NA
# temp = data.frame(SampleID = taxa$SampleID, Response = resp)
# taxa = merge(temp, taxa, by='SampleID')

# ## write dataframe to file
write.table( taxa, file = out_file, sep = '\t',quote = F,row.names = F)


