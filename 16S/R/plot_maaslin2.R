source('/pita/users/tzipi/code/R_figs/maaslin2_funs.R')
source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
library(ggplot2)

fdr = 0.25
taxonomy_flag = T
taxonomy_file = '/pita/users/tzipi/projects/16S/all_merges/DB1-20_merge/res/16S_DB1-20_merged/taxonomy/metadata_v2_source.tsv'

name = 'source'
path = '../res_china_stool_amnonFlt/'
# path = '../res_israel_stool_amnonFlt/'
# path = '../res_israel_stool_amnonFlt/'
# path = '../res_israel_noRS_amnonFlt/'
taxa_file = sprintf('%s/biom/%s_deblur_map_SVs.1.txt',path, name)
na_str = c('no_data','_','NA','unknown', 'other','na','No_followup')
taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str)
head_path = sprintf('%s/R_res/maaslin2/',path)

# name = 'israel_noRS_CD_age_gender_location2_pnIDrnd'
# taxa$location[taxa$location %in% c('R','TI')] = 'Biopsy'
# taxa$Patient_group = sprintf('%s %s', taxa$Patient_group2, taxa$location)
# x_var = 'Age_years'; group_var = 'Patient_group'; taxa_sig_var = 'Patient_group2'

# name = 'urban_CD_age_gender'
# taxa = taxa[taxa$Patient_group %in% c('Chinese_crohns','Urban_health'),]
# x_var = 'Age_years'; group_var = 'Patient_group'; taxa_sig_var = 'Patient_group2'

name = 'controls_Fat_gr_age_gender_group2RuralUnderBase'
x_var = 'Fat_gr_day'; group_var = 'Patient_group2'; taxa_sig_var = 'Fat_gr_day'

# name = 'controls_Ferrus_age_gender_group2RuralUnderBase'
# x_var = 'Ferrus_mg_day'; group_var = 'Patient_group2'; taxa_sig_var = 'Ferrus_mg_day'

taxa_m = maaslin2_set_data_for_heatmap_plot(taxa, head_path, name, x_var, 
                                            taxa_sig_var, group_var = group_var, 
                                            taxonomy_file = taxonomy_file, fdr = fdr,
                                            taxonomy_flag = taxonomy_flag)
library(viridis)
taxa_m$Dx = ifelse(grepl(pattern = '[Cc]rohns',x = taxa_m$Patient_group),'CD','Control')
# taxa_m$dir = ifelse(taxa_m$coef > 0,'Down','Up')
taxa_m$dir = ifelse(taxa_m$coef > 0,'Up','Down')
taxa_m$dir = factor(taxa_m$dir, levels= c('Up','Down'))

taxa_m$Group = taxa_m$Dx; 
taxa_m$Group[taxa_m$Patient_group2 == 'Rural_health_<50%_in_city'] = 'Rural'
taxa_m$Group[taxa_m$Patient_group2 == 'Rural_health_>50%_in_city']= 'Rural-\nUrban'
taxa_m$Group[taxa_m$Patient_group2 == 'Urban_health']= 'Urban'
# taxa_m$Group = factor(taxa_m$Group, levels = c('CD','Urban','Rural\n>50%\nin city','Rural\n<50%\nin city'))
taxa_m$Group = factor(taxa_m$Group, levels = c('CD','Urban','Rural-\nUrban','Rural'))
p_heatmap = ggplot(taxa_m, aes_string(x=x_var, y='taxa', fill = 'log10_RA')) + geom_tile() + 
  theme_minimal() + 
  ylab('ASV') + xlab('') + 
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.spacing = unit(0.1,'lines'), plot.background = element_rect(fill = "white"),
        text = element_text(size=16), strip.text.x = element_text(angle = 45)) +
  scale_fill_viridis() +
  facet_grid(dir~Group, scales = 'free', space = 'free', switch = "both") 
  # facet_grid(dir~Dx, scales = 'free', space = 'free', switch = "both") 
  # facet_grid(dir~Dx+loc, scales = 'free', space = 'free', switch = "both") 

out_path = 'art_plots/res/'
out_fig = sprintf('%s/%s_maaslin2_heatmap_fdr%s.tiff', out_path, name, fdr)
# ggsave(out_fig,plot = p_heatmap, device = 'tiff', width = 5,height = 3, compression ='lzw', dpi=600)
ggsave(out_fig,plot = p_heatmap, device = 'tiff', width = 5,height = 3.5, compression ='lzw', dpi=600)

# path = '../res_israel_noRS_amnonFlt/'
# taxa_file = sprintf('%s/biom/%s_deblur_map_SVs.1.txt',path, 'source')
# taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str)
# # p_heatmap25 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', add_name = '_israel_noRS', taxa_sig_var = taxa_sig_var, fdr = 0.25)
# # p_heatmap1 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', add_name = '_israel_noRS', taxa_sig_var = taxa_sig_var, fdr = 0.1)
