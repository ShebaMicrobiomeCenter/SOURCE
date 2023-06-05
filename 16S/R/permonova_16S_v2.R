library(ggplot2)
library(vegan)
source('/pita/users/tzipi/code/R_figs/figs_funcs.R')

source('/pita/users/tzipi/code/R_figs/permanova_funcs.R')

set.seed(4)
# set.seed(5)

## uses permanova to calculate the microbiome varience explained by metadata.
## R2 is the varience explained in the adonis2 result. 

country = 'china'
country = 'israel'

main_name = 'source_16S_un_v2'
# main_name = sprintf('16S_un_%s', country)

type = 'unweighted_unifrac'

age_match_flag = F

path = sprintf('../../16S/res_%s_stool_amnonFlt/', country)
# path = '../../16S/res_israel_stool_amnonFlt/'

# map_file = sprintf('%s/map.txt', path)
map_file = '../../metadata/alona_data/SOURCE_Israel_China_data_v11.txt'
distance_matrix_file = sprintf('%s/core-metrics-results/%s_distance_matrix/distance-matrix.tsv',path, type)
out_path =  sprintf('res/')

# dir.create(out_path)    

extra_vars = c('Age_years','Gender')
# extra_vars = c()
strata_var = NA

na_str = c('no_data','_','NA','unknown', 'other','na','No_followup','ND','Unkown')
dist = read.table(distance_matrix_file, header = T, row.names = 1, sep = '\t')
map = read.table(map_file, header = T, sep = '\t', row.names = 1, na.strings = na_str)

if ( age_match_flag )
{
  map = map[map$Age_years >= 25 & map$Age_years <= 53,]
  main_name = sprintf('%s_ageMatch', main_name)
}

# ffq_data = '../../metadata/FFQ/FFQ_v9_kcal_norm.txt'
# ffq_data = '../../metadata/FFQ/FFQ_v9_absolute.txt'
ffq_data = '../../metadata/FFQ/FFQ_v10_absolute.txt'
ffq = read.table(ffq_data)
ffq$pn_ID = row.names(ffq)

map = map[,c('pn_ID',names(map)[!names(map) %in% names(ffq) ])]
map$sampleID = row.names(map)
map2 =merge(map, ffq, by='pn_ID',all.x = T)
row.names(map2) = map2$sampleID
map = map2

map_all = map
dist_all = dist

# name = 'China_stool_all'
# map = map[ map$Cohort == 'SOURCE_China',]
# map = map[ map$location == 'Stool',]

name = 'China_stool_CD'
map = map[ map$Cohort == 'SOURCE_China',]
map = map[ map$Patient_group == 'Chinese_crohns',]
map = map[ map$location == 'Stool',]

# name = 'China_stool_Healthy'
# map = map[ map$Cohort == 'SOURCE_China',]
# map = map[ map$Patient_group %in% c('Urban_health','Rural_health') , ]
# map = map[ map$location == 'Stool',]

# name = 'China_stool_Healthy_Urban'
# map = map[ map$Cohort == 'SOURCE_China',]
# map = map[ map$Patient_group %in% c('Urban_health') , ]
# map = map[ map$location == 'Stool',]

# name = 'China_stool_Healthy_Rural'
# map = map[ map$Cohort == 'SOURCE_China',]
# map = map[ map$Patient_group %in% c('Rural_health') , ]
# map = map[ map$location == 'Stool',]

# name = 'China_stool_Healthy_Rural>50%'
# map = map[ map$Cohort == 'SOURCE_China',]
# map = map[ map$Patient_group2 %in% c('Rural_health_>50%_in_city') , ]
# map = map[ map$location == 'Stool',]

# name = 'China_stool_Healthy_Rural<50%'
# map = map[ map$Cohort == 'SOURCE_China',]
# map = map[ map$Patient_group2 %in% c('Rural_health_<50%_in_city') , ]
# map = map[ map$location == 'Stool',]

# name = 'Israel_stool_all'
# map = map[ map$Cohort == 'SOURCE_Israel',]
# map = map[ map$location == 'Stool',]

# name = 'Israel_stool_CD'
# map = map[ map$Cohort == 'SOURCE_Israel',]
# map = map[ map$Patient_group == 'Israeli_Crohns',]
# map = map[ map$location == 'Stool',]

# name = 'Israel_stool_Healthy'
# map = map[ map$Cohort == 'SOURCE_Israel',]
# map = map[ map$Patient_group %in% c('Israeli_healthy') , ]
# map = map[ map$location == 'Stool',]

vars_2_check = names(map)
perm_df = permanova_organize_and_run( map, dist, vars_2_check, min_notNA_num = 3, extra_vars, strata_var)
perm_df$Type = name

out_path = sprintf('%s/%s/',out_path, main_name)
dir.create(out_path)
out_file = sprintf('%s/%s_%s.txt', out_path, name, perm_df$Permanova_extra_prms[1])
write.table(x=perm_df, file = out_file, quote = F,sep = '\t',row.names = F)



