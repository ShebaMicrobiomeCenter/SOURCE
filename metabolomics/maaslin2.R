source('/pita/users/tzipi/code/R_figs/maaslin2_funs.R')
source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
library(Maaslin2)
library(ggplot2)

taxonomy_flag = F
# taxonomy_file = '/pita/users/tzipi/projects/16S/all_merges/DB1-20_merge/res/16S_DB1-20_merged/taxonomy/metadata_v2_source.tsv'

name = 'source'
path = '../data/'
# ftr_file = sprintf('%s/SOURCE2022_%s_Data_MZM13_N2_filtered_calour.csv',path, type)
ftr_file = '../data/SOURCE_china_mets_table_norm_v2.tsv'
map_file = '../data/SOURCE_china_mets_map.tsv'

na_str = c('no_data','_','NA','unknown', 'other','na','No_followup')
ftr = read.table(file = ftr_file, header = T,row.names = 1, na.strings = na_str, sep='\t')
map = read.table(map_file,sep="\t", header=TRUE, na.strings = na_str)

ftr = as.data.frame(t(ftr))

head_path = sprintf('res/maaslin2_v2/')
dir.create(head_path)

map$Description = map$SampleID
map$SampleID = make.names(map$SampleID)
ftr$SampleID = make.names(row.names(ftr))
taxa = merge(map, ftr, by = 'SampleID')

# name = sprintf('Dx_gender_age')
# fixed_cols = c('Dx', 'Gender' ,'Age_years')
# # fixed_cols = c('Dx', 'Gender' ,'Age')
# # random_cols = c('pn_ID')
# random_cols = c()
# metadata_cols = c(fixed_cols, random_cols)
# x_var = 'Age_years'; group_var = 'Dx'; taxa_sig_var = 'Dx'

# name = sprintf('Dx_gender_age_group')
# fixed_cols = c('Dx', 'Gender' ,'Age_years','Patient_group2')
# # fixed_cols = c('Dx', 'Gender' ,'Age')
# # random_cols = c('pn_ID')
# random_cols = c()
# levels_ref = c('Patient_group2,Rural_health_<50%_in_city')
# metadata_cols = c(fixed_cols, random_cols)
# x_var = 'Age_years'; group_var = 'Dx'; taxa_sig_var = 'Dx'

# name = sprintf('gender_age_group')
# fixed_cols = c('Gender' ,'Age_years','Patient_group2')
# # fixed_cols = c('Dx', 'Gender' ,'Age')
# # random_cols = c('pn_ID')
# random_cols = c()
# levels_ref = c('Patient_group2,Rural_health_<50%_in_city')
# metadata_cols = c(fixed_cols, random_cols)
# x_var = 'Age_years'; group_var = 'Dx'; taxa_sig_var = 'Dx'

# taxa = taxa[taxa$Patient_group2 %in% c('Urban_health','Chinese_crohns'),]
# name = sprintf('gender_age_group_inUrban')
# fixed_cols = c('Gender' ,'Age_years','Patient_group2')
# # fixed_cols = c('Dx', 'Gender' ,'Age')
# # random_cols = c('pn_ID')
# random_cols = c()
# levels_ref = c('Patient_group2,Rural_health_<50%_in_city')
# metadata_cols = c(fixed_cols, random_cols)
# x_var = 'Age_years'; group_var = 'Dx'; taxa_sig_var = 'Dx'

taxa = taxa[! taxa$Patient_group2 %in% c('Urban_health','Chinese_crohns'),]
name = sprintf('gender_age_group_inRural')
fixed_cols = c('Gender' ,'Age_years','Patient_group2')
# fixed_cols = c('Dx', 'Gender' ,'Age')
# random_cols = c('pn_ID')
random_cols = c()
levels_ref = c('Patient_group2,Rural_health_<50%_in_city')
metadata_cols = c(fixed_cols, random_cols)
x_var = 'Age_years'; group_var = 'Dx'; taxa_sig_var = 'Dx'


maas_path = sprintf('%s/%s',head_path, name)

ASV = write_maaslin2(taxa, metadata_cols, maas_path, taxonomy_flag = taxonomy_flag, taxonomy_file = taxonomy_file)
# ASV = write_maaslin2(taxa, metadata_cols, maas_path, taxonomy_flag = taxonomy_flag)

out_path = sprintf('%s/res/',maas_path)

fit_data <- Maaslin2( input_data = sprintf('%s/taxonomy.tsv',maas_path),
                      input_metadata = sprintf('%s/metadata.tsv',maas_path),
                      output = out_path,
                      reference = levels_ref,
                      # transform = "AST",
                      # random_effects = random_cols,
                      fixed_effects = fixed_cols)


# # if changes to taxonomy, add the asv to the result files
# if ( taxonomy_flag )
# {
#   maaslin_add_asv_to_res(ASV, out_path)
# }
# maaslin2_plot_colored_scatters(taxa, x_var = taxa_sig_var, col_var = group_var, head_path, name, out_dir_name = sprintf('my_%s_scatters', taxa_sig_var), fdr = 0.25)

# p_heatmap25 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', taxa_sig_var = taxa_sig_var, fdr = 0.25, taxonomy_flag = F)
# p_heatmap1 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', taxa_sig_var = taxa_sig_var, fdr = 0.1)
# # taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str)
# # p_heatmap25 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', add_name = '_all_samples', taxa_sig_var = taxa_sig_var, fdr = 0.25)
# # p_heatmap1 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', add_name = '_all_samples', taxa_sig_var = taxa_sig_var, fdr = 0.1)
# # 

# path = '../res_israel_noRS_amnonFlt/'
# taxa_file = sprintf('%s/biom/%s_deblur_map_SVs.1.txt',path, 'source')
# taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str)
# # p_heatmap25 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', add_name = '_israel_noRS', taxa_sig_var = taxa_sig_var, fdr = 0.25)
# # p_heatmap1 = maaslin2_healtmap_plot(taxa, head_path, name, x_var, group_var = group_var, taxonomy_file, out_dir_name = 'my_heatmaps', add_name = '_israel_noRS', taxa_sig_var = taxa_sig_var, fdr = 0.1)

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

