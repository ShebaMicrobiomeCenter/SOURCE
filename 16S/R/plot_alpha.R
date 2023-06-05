library(ggplot2)
# library(stringr)

source('/pita/users/tzipi/code/R_figs/figs_funcs.R')

test_type = 'wilcox'

name = 'source'

# path = '../res_israel_stool_amnonFlt/'
path = '../res_israel_noRS_amnonFlt/'
# path = '../res_china_stool_amnonFlt/'
# path = '../res_china_stool_rural_ageMatched_amnonFlt/'
taxa_file = sprintf('%s/biom/%s_deblur_map_L7.1.txt',path, name)
out_path =  sprintf('%s/R_res',path)
dir.create(out_path)

na_str = c('no_data','_','NA','unknown', 'other','na','No_followup','ND')
taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str)


fdr = 'BH' # 'bonferroni'
alpha_var = 'faith_pd'; y_lb = 'Faith\'s Phylogenetic Diversity'

check_val = 'Patient_group2'

# taxa$Patient_group[taxa$Patient_group == 'Chinese_crohns'] = 'A Chinese_crohns'
# taxa$Patient_group[taxa$Patient_group == 'Urban_health'] = 'B Urban_health'
# taxa$Patient_group[taxa$Patient_group == 'Israeli_Crohns'] = 'D Israeli_Crohns'
# taxa$Patient_group[taxa$Patient_group == 'Israeli_healthy'] = 'E Israeli_healthy'
# taxa$Patient_group2[taxa$Patient_group == 'Rural_health'] = 'C Rural_health'

taxa$Patient_group2[taxa$Patient_group2 == 'Chinese_crohns'] = 'China CD'
taxa$Patient_group2[taxa$Patient_group2 == 'Urban_health'] = 'China Urban'
taxa$Patient_group2[taxa$Patient_group2 == 'Rural_health_>50%_in_city'] = 'Rural-Urban'
taxa$Patient_group2[taxa$Patient_group2 == 'Rural_health_<50%_in_city'] = 'Rural'

taxa$Patient_group[taxa$Patient_group == 'Chinese_crohns'] = 'A Chinese_crohns'
taxa$Patient_group[taxa$Patient_group == 'Urban_health'] = 'B Urban_health'
taxa$Patient_group[taxa$Patient_group == 'Rural_health'] = 'C Rural_health'

taxa2 = taxa

name = 'isreal_noRS'

# name = 'isreal_stool'
# taxa2 = taxa2[taxa2$location == 'Stool',]

# name = 'isreal_TI'
# taxa2 = taxa2[taxa2$location == 'TI',]

# name = 'isreal_R'
# taxa2 = taxa2[taxa2$location == 'R',]

# name = 'china'

# name = 'china_rural'
# taxa2 = taxa2[taxa2$Patient_group == 'C Rural_health',]

# name = 'china_urban_CD'
# taxa2 = taxa2[taxa2$Patient_group %in% c('A Chinese_crohns','B Urban_health'),]

out_path = 'art_plots/res/'
alpha_res2 = plot_variable_and_stats(taxa2, check_val, y_val =alpha_var, y_lb = y_lb, sig_ast_flag = T, print_pvals = T, boxplot_p_text_flag = F, jitter_flag = F, boxplot_test_type = 'wilcox')
alpha_res2 = alpha_res2 + 
  scale_fill_manual(values= c('gray80','gray80','gray80','gray80')) + 
  xlab('') + theme(panel.grid = element_blank())
ggsave(sprintf( '%s/%s_%s.tiff', out_path, sprintf('Alpha_Diversity_%s', check_val), name ),
       plot = alpha_res2, device = 'tiff', width =2,height = 3.5, compression = 'lzw',dpi=600)



