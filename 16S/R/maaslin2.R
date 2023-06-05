source('/pita/users/tzipi/code/R_figs/maaslin2_funs.R')
source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
library(Maaslin2)
library(ggplot2)

taxonomy_flag = T
taxonomy_file = '/pita/users/tzipi/projects/16S/all_merges/DB1-20_merge/res/16S_DB1-20_merged/taxonomy/metadata_v2_source.tsv'

prms = c('SampleID','dysbiosis_index')

name = 'source'
# path = '../res_china_stool_amnonFlt/'
# path = '../res_israel_stool_amnonFlt/'
# path = '../res_israel_stool_amnonFlt/'
path = '../res_israel_noRS_amnonFlt/'
taxa_file = sprintf('%s/biom/%s_deblur_map_SVs.1.txt',path, name)
na_str = c('no_data','_','NA','unknown', 'other','na','No_followup')
taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str)

head_path = sprintf('%s/R_res/maaslin2/',path)
dir.create(head_path)

# name = 'controls_Ferrus_age_gender_group2RuralUnderBase'
# taxa = taxa[taxa$Dx == 'healthy',]
# fixed_cols = c('Ferrus_mg_day','Age_years','Gender','Patient_group2')
# # levels_ref = c("Status,control;any_treatment_before_YN,control;Type_Treatment,control;Histology,control;Stage_at_treatment_start,control")
# levels_ref = c("Patient_group2,Rural_health_<50%_in_city")
# # random_cols = c('Cage.x')
# random_cols = c()
# metadata_cols = c(fixed_cols, random_cols)
# x_var = 'Ferrus_mg_day'; group_var = 'Patient_group2'; taxa_sig_var = 'Ferrus_mg_day'

# name = 'controls_Fat_gr_age_gender_group2RuralUnderBase'
# taxa = taxa[taxa$Dx == 'healthy',]
# fixed_cols = c('Fat_gr_day','Age_years','Gender','Patient_group2')
# # levels_ref = c("Status,control;any_treatment_before_YN,control;Type_Treatment,control;Histology,control;Stage_at_treatment_start,control")
# levels_ref = c("Patient_group2,Rural_health_<50%_in_city")
# # random_cols = c('Cage.x')
# random_cols = c()
# metadata_cols = c(fixed_cols, random_cols)
# x_var = 'Fat_gr_day'; group_var = 'Patient_group2'; taxa_sig_var = 'Fat_gr_day'

# name = 'controls_Saturated_fat_gr_age_gender_group2RuralUnderBase'
# taxa = taxa[taxa$Dx == 'healthy',]
# fixed_cols = c('Saturated_fat_gr_day','Age_years','Gender','Patient_group2')
# # levels_ref = c("Status,control;any_treatment_before_YN,control;Type_Treatment,control;Histology,control;Stage_at_treatment_start,control")
# levels_ref = c("Patient_group2,Rural_health_<50%_in_city")
# # random_cols = c('Cage.x')
# random_cols = c()
# metadata_cols = c(fixed_cols, random_cols)
# x_var = 'Saturated_fat_gr_day'; group_var = 'Patient_group2'; taxa_sig_var = 'Saturated_fat_gr_day'

# name = 'rural_age_gender'
# taxa = taxa[taxa$Patient_group == 'Rural_health',]
# fixed_cols = c('Age_years','Gender','Patient_group2')
# # levels_ref = c("Status,control;any_treatment_before_YN,control;Type_Treatment,control;Histology,control;Stage_at_treatment_start,control")
# # levels_ref = c("Patient_group2,Rural_health_<50%_in_city")
# # random_cols = c('Cage.x')
# random_cols = c()
# metadata_cols = c(fixed_cols, random_cols)
# x_var = 'Age_years'; group_var = 'Patient_group2'; taxa_sig_var = 'Patient_group2'

# name = 'urban_CD_age_gender'
# taxa = taxa[taxa$Patient_group %in% c('Chinese_crohns','Urban_health'),]
# fixed_cols = c('Age_years','Gender','Patient_group2')
# # levels_ref = c("Status,control;any_treatment_before_YN,control;Type_Treatment,control;Histology,control;Stage_at_treatment_start,control")
# # levels_ref = c("Patient_group2,Rural_health_<50%_in_city")
# # random_cols = c('Cage.x')
# random_cols = c()
# metadata_cols = c(fixed_cols, random_cols)
# x_var = 'Age_years'; group_var = 'Patient_group2'; taxa_sig_var = 'Patient_group2'

# name = 'israel_stool_CD_age_gender'
# fixed_cols = c('Age_years','Gender','Patient_group2')
# # levels_ref = c("Status,control;any_treatment_before_YN,control;Type_Treatment,control;Histology,control;Stage_at_treatment_start,control")
# # levels_ref = c("Patient_group2,Rural_health_<50%_in_city")
# # random_cols = c('Cage.x')
# random_cols = c()
# metadata_cols = c(fixed_cols, random_cols)
# x_var = 'Age_years'; group_var = 'Patient_group2'; taxa_sig_var = 'Patient_group2'

name = 'israel_noRS_CD_age_gender_location2_pnIDrnd'
fixed_cols = c('Age_years', 'Gender', 'Patient_group2', 'location')
taxa$location[taxa$location %in% c('R','TI')] = 'Biopsy'
# levels_ref = c("Status,control;any_treatment_before_YN,control;Type_Treatment,control;Histology,control;Stage_at_treatment_start,control")
# levels_ref = c("Patient_group2,Rural_health_<50%_in_city")
# random_cols = c('Cage.x')
random_cols = c('pn_ID')
metadata_cols = c(fixed_cols, random_cols)
taxa$Patient_group = sprintf('%s %s', taxa$Patient_group2, taxa$location)
x_var = 'Age_years'; group_var = 'Patient_group'; taxa_sig_var = 'Patient_group2'

maas_path = sprintf('%s/%s',head_path, name)

ASV = write_maaslin2(taxa, metadata_cols, maas_path, taxonomy_flag = taxonomy_flag, taxonomy_file = taxonomy_file)
# ASV = write_maaslin2(taxa, metadata_cols, maas_path, taxonomy_flag = taxonomy_flag)

out_path = sprintf('%s/res/',maas_path)

fit_data <- Maaslin2( input_data = sprintf('%s/taxonomy.tsv',maas_path),
                      input_metadata = sprintf('%s/metadata.tsv',maas_path),
                      output = out_path,
                      # reference = levels_ref,
                      # transform = "AST",
                      fixed_effects = fixed_cols,
                      random_effects = random_cols)

# if changes to taxonomy, add the asv to the result files
if ( taxonomy_flag )
{
  maaslin_add_asv_to_res(ASV, out_path)
}
# maaslin2_plot_colored_scatters(taxa, x_var = taxa_sig_var, col_var = group_var, head_path, name, out_dir_name = sprintf('my_%s_scatters', taxa_sig_var), fdr = 0.25)

p_heatmap25 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', taxa_sig_var = taxa_sig_var, fdr = 0.25)
p_heatmap1 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', taxa_sig_var = taxa_sig_var, fdr = 0.1)
# taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str)
# p_heatmap25 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', add_name = '_all_samples', taxa_sig_var = taxa_sig_var, fdr = 0.25)
# p_heatmap1 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', add_name = '_all_samples', taxa_sig_var = taxa_sig_var, fdr = 0.1)
# 

path = '../res_israel_noRS_amnonFlt/'
taxa_file = sprintf('%s/biom/%s_deblur_map_SVs.1.txt',path, 'source')
taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str)
# p_heatmap25 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', add_name = '_israel_noRS', taxa_sig_var = taxa_sig_var, fdr = 0.25)
# p_heatmap1 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', add_name = '_israel_noRS', taxa_sig_var = taxa_sig_var, fdr = 0.1)

# fdr = 0.25
# taxa_m = maaslin2_set_data_for_heatmap_plot(taxa, head_path, name, x_var, taxa_sig_var, group_var = group_var, taxonomy_file = taxonomy_file, fdr)
# taxa_m = merge(taxa_m, taxa[,c('SampleID','location')], by = 'SampleID', all.x=T, sort = F)
# 
# library(viridis)
# p_heatmap = ggplot(taxa_m, aes_string(x=x_var, y='taxa', fill = 'log10_RA')) + geom_tile() + 
#   facet_grid(~Patient_group2+location, scales = 'free') + 
#   theme_classic() + 
#   ylab('ASV') + 
#   scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#   scale_fill_viridis() 
# out_path = sprintf('%s/%s/res/%s',head_path,name, 'my_heatmaps')
# dir.create(out_path)
# out_fig = sprintf('%s/%s_sigTaxa_%s_%s_heatmap%s_fdr%s.tiff', out_path, taxa_sig_var, x_var, group_var, '_israel_noRS', fdr)
# ggsave(out_fig,plot = p_heatmap, device = 'tiff', width = 12,height = 2+length(unique(taxa_m$ASV))/8, compression  = 'lzw')

