source('/pita/users/tzipi/code/R_figs/maaslin2_funs.R')
source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
library(Maaslin2)
library(ggplot2)
library(patchwork)

taxonomy_file = '/pita/users/tzipi/projects/16S/all_merges/DB1-20_merge/res/16S_DB1-20_merged/taxonomy/metadata_v2_source.tsv'

name = 'source'
# path = '../res_china_stool_amnonFlt/'
path = '../res_israel_noRS_amnonFlt/'
taxa_file = sprintf('%s/biom/%s_deblur_map_SVs.1.txt',path, name)
bdiv_file = sprintf('%s/core-metrics-results/unweighted_unifrac_distance_matrix/distance-matrix.tsv',path)
na_str = c('no_data','_','NA','unknown', 'other','na')
taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str)
taxa$location[taxa$location %in% c('R','TI')] = 'Biopsy'
taxa$sample_ID = gsub('-','_',taxa$sample_ID)

bdiv = read.table(bdiv_file, header = T, row.names = 1, sep = '\t')

pos = !grepl('_',taxa$sample_ID)
taxa$sample_ID[pos] = sprintf('%s_%s', taxa$pn_ID[pos], taxa$Location_inflamation[pos])

# add_name = ''

# add_name = '_healthy'
# taxa = taxa[taxa$Dx == 'healthy',]

add_name = '_Crohns'
taxa = taxa[taxa$Dx == 'Crohns',]


sv_df = taxa[,(which(names(taxa)=='Description')+1):dim(taxa)[2] ]
row.names(sv_df) = taxa$sample_ID
sv_df = as.data.frame(t(sv_df))

# sv_df = taxa[,c(1,(which(names(taxa)=='Description')+1):dim(taxa)[2]) ]
# sv_df = reshape2::melt(sv_df, value.name = 'RA')
# names(sv_df)[2] = 'SV'
# 
# 
# sv_df = merge(sv_df, taxa[,c('SampleID','pn_ID','Cohort','Patient_group','Dx','Location_inflamation','location')])

dist_type = 'bDiv'; dist_lab = 'Unweighted Unifrac distance\nbetween 2 samples'
# dist_type = 'corr'; dist_lab = 'Spearman correlations r\nof all taxa between 2 samples'

s1 = c()
s2 = c()
r = c()
p = c()
smps = unique(taxa$sample_ID)
for ( i in 1:length(smps) )
{
  for (j in i:length(smps))
  {
    if ( i != j )
    {
      # cor_res = cor.test(x=sv_df$RA[sv_df$SampleID == smps[i]], y= sv_df$RA[sv_df$SampleID == smps[j]], method = 'spearman')
      cor_res = cor.test(sv_df[, smps[i]],sv_df[, smps[j]], method = 'spearman')
      # g = ggplot(sv_df, aes_string(x=smps[i], y=smps[j])) + geom_point() + 
      #   scale_x_log10() + scale_y_log10() + theme_bw() +
      #   geom_abline(slope=1, intercept=0, colour = 'gray') + 
      #   ggtitle(sprintf('Spearman r = %.2f, p = %.2e', cor_res$estimate, cor_res$p.value))
      s1 = c(s1, smps[i]) 
      s2 = c(s2, smps[j]) 
      if (dist_type == 'bDiv')
      {
        r = c(r, bdiv[ taxa$SampleID[taxa$sample_ID == smps[i]], taxa$SampleID[taxa$sample_ID == smps[j]] ])
      } else if (dist_type == 'corr')
        r = c(r, cor_res$estimate)
      # p = c(p, cor_res$p.value)
    }
  }
}
cor_df = data.frame(x=s1, y=s2, r=r)

cor_df$x_pnID = gsub('_.*','',cor_df$x)
cor_df$y_pnID = gsub('_.*','',cor_df$y)
cor_df$x_location_inflamation = gsub('.*_','',cor_df$x)
cor_df$y_location_inflamation = gsub('.*_','',cor_df$y)
cor_df$x_location_type = gsub('[1-4]','',cor_df$x_location_inflamation)
cor_df$y_location_type = gsub('[1-4]','',cor_df$y_location_inflamation)
cor_df$x_location = ifelse(cor_df$x_location_inflamation == 'Stool','Stool','Biopsy')
cor_df$y_location = ifelse(cor_df$y_location_inflamation == 'Stool','Stool','Biopsy')

cor_df$locations = sprintf('%s x %s', cor_df$x_location_type, cor_df$y_location_type)
cor_df$locations_same = ifelse(cor_df$x_location == cor_df$y_location,'Same location','different locations')
cor_df$Same_patient = ifelse(cor_df$x_pnID==cor_df$y_pnID, 'Same patient','Different patients')
cor_df$same_pnID = as.factor(cor_df$x_pnID==cor_df$y_pnID)
cor_df$same_location = as.factor(cor_df$x_location == cor_df$y_location)

cor_df$locations[cor_df$locations %in% c('R x Stool','Stool x R')] = 'Stool x R'
cor_df$locations[cor_df$locations %in% c('Stool x TI','TI x Stool')] = 'Stool x TI'
cor_df$locations[cor_df$locations %in% c('TI x R','R x TI')] = 'TI x R'

cor_df2 = cor_df
cor_df2 = cor_df2[cor_df2$locations %in% c('Stool x R','Stool x TI','TI x R'),]
cor_df2$locations2 = ifelse(cor_df2$locations %in% c('TI x R'),
                            'Different biopsies\n(R x TI)','Stool x Biopsy')
cor_df2$locations2 = factor(cor_df2$locations2, levels = c('Stool x Biopsy','Different biopsies\n(R x TI)'))
g = ggplot(cor_df2) + 
  geom_boxplot(aes(y=r, x=Same_patient, fill = locations2)) + 
  theme_classic() + ylab(dist_lab) + 
  # xlab('Samples from the same patient')+ 
  xlab('') + 
  scale_fill_manual(values = c('gray50','gray80'), name = '')
res_g = boxplot_2_vars_grouped_add_asterix( g, cor_df2, 
      x_Val ='Same_patient', y_val = 'r', fill_val ='locations2', 
      test_type = 'wilcox', fdr= 3)
out_path = 'res/location_correlations/'
dir.create(out_path)
out_fig = sprintf('%s/%s_biopsy_stool_dist_boxplot.tiff',out_path, dist_type)
res_g = format_fig(res_g, bw_flag = F, set_ttl_flag = F, x_angle = 45)
ggsave(out_fig,plot = res_g, device = 'tiff', width = 6,height = 5, 
       compression  = 'lzw', dpi = 300)


g2 = ggplot(cor_df2, aes(x=locations, y=r)) + 
  geom_boxplot(fill =  'gray') + 
  facet_wrap(~Same_patient, scales = 'free_x') + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid = element_blank()) + 
  ggsignif::geom_signif( comparisons = list(c("Stool x R", "Stool x TI"),
                                            c("Stool x R", "TI x R"),
                                            c("Stool x TI", "TI x R")),
                         step_increase = 0.1, map_signif_level = T, test = 'wilcox.test') + 
  coord_cartesian(ylim = c(0.2,1)) + xlab('') + 
  ylab(dist_lab) 


# t = 'Same patient'
# gp1 = plot_variable_and_stats(taxa = cor_df2[cor_df2$Same_patient == t,],
#                               check_val = 'locations', y_val = 'r', 
#                               sig_ast_flag = T,print_pvals = F, y_lb = '',
#                               jitter_flag = F) + ggtitle(t)
# # gp1 = ggplot(cor_df2[cor_df2$Same_patient == t,], aes(x=locations, y=r)) + 
# #   geom_boxplot(fill =  'gray') + 
# #   facet_wrap(~Same_patient, scales = 'free_x') + 
# #   theme_bw() + 
# #   theme(axis.text.x = element_text(angle = 45, hjust = 1), 
# #         panel.grid = element_blank()) + 
# #   ggsignif::geom_signif( comparisons = list(c("Stool x R", "Stool x TI"),
# #                                             c("Stool x R", "TI x R"),
# #                                             c("Stool x TI", "TI x R")),
# #                          step_increase = 0.1, map_signif_level = T, test = 'wilcox.test') + 
# #   coord_cartesian(ylim = c(0.2,1)) + xlab('') + 
# #   ylab('Unweighted Unifrac distance\nbetween 2 samples') 
# gp1 = ggplot(cor_df2[cor_df2$Same_patient == t,]) + 
#   geom_boxplot(fill =  'gray', aes(x=locations, y=r)) + 
#   # facet_wrap(~Same_patient, scales = 'free_x') + 
#   ggtitle(t) + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), 
#         panel.grid = element_blank()) + 
#   xlab('') + 
#   ylab(dist_lab) 
# res = add_significant_asterix_to_plot_v3(p = gp1, 
#       var = as.factor(cor_df2$locations[cor_df2$Same_patient == t]),
#       val = cor_df2$r[cor_df2$Same_patient == t],
#       test_type = 'wilcox', print_pvals = T,  fdr_flag = F
#         )
# gp1 = res[[1]] + ggtitle(t)
# 
# t2 = 'Different patients'
# # gp2 = plot_variable_and_stats(taxa = cor_df2[cor_df2$Same_patient == t2,],
# #                               check_val = 'locations', y_val = 'r', 
# #                               sig_ast_flag = T,print_pvals = F, y_lb = '',
# #                               jitter_flag = F, boxplot_test_type = 'wilcox')+ ggtitle(t2)
# 
# gp2 = ggplot(cor_df2[cor_df2$Same_patient == t2,]) + 
#   geom_boxplot(fill =  'gray', aes(x=locations, y=r)) + 
#   # facet_wrap(~Same_patient, scales = 'free_x') + 
#   ggtitle(t) + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), 
#         panel.grid = element_blank()) + 
#   xlab('') + 
#   ylab(dist_lab) 
# res = add_significant_asterix_to_plot_v3(p = gp2, 
#                                          var = as.factor(cor_df2$locations[cor_df2$Same_patient == t2]),
#                                          val = cor_df2$r[cor_df2$Same_patient == t2],
#                                          test_type = 'wilcox', print_pvals = T,  fdr_flag = F
# )
# gp2 = res[[1]] + ggtitle(t2)
# 
# gp = (gp2 + ylim(c(0.2,0.9))) + (gp1 + ylab('') + ylim(c(0.2,0.9)))
# gp
# 
# cor_df2 = cor_df
# cor_df2 = cor_df2[cor_df2$x_location !='Stool' & cor_df2$y_location !='Stool',]
# cor_df2$xy_location_inflamation = 'temp'
# for ( i in 1:length(cor_df2$x_location_inflamation))
# {
#   temp = c(cor_df2$x_location_inflamation[i],cor_df2$y_location_inflamation[i])
#   temp = sort(temp)
#   cor_df2$xy_location_inflamation[i] = sprintf('%s_%s',temp[1],temp[2])
# }
# 
# 
# cor_df3 = cor_df2
# cor_df3 = cor_df3[cor_df3$same_pnID == 'TRUE',]
# res2 = plot_variable_and_stats(taxa = cor_df3, check_val = 'xy_location_inflamation',y_val = 'r',sig_ast_flag = T, print_pvals = T, y_lb = 'Unweighted Unifrac distance\nbetween 2 samples')
# out_fig = sprintf('%s/2_patients_taxa_bdiv_withtinBiopsy_samePatient_location_inflamation%s.tiff',out_path, add_name)
# # ggsave(out_fig,plot = res2, device = 'tiff', width = 6,height = 4, compression  = 'lzw', dpi = 300)
# 
# 
# cor_df4 = cor_df2
# cor_df4 = cor_df4[cor_df4$xy_location_inflamation %in% c('R3_R4','TI1_TI2','R4_TI2'),]
# cor_df4$xy_location_inflamation = factor(cor_df4$xy_location_inflamation,levels = c('R3_R4','TI1_TI2','R4_TI2'))
# g4 = ggplot(cor_df4) + geom_boxplot(aes(y=r, x=same_pnID , fill = xy_location_inflamation)) + 
#   theme_bw()  +  ylab('Unweighted Unifrac distance\nbetween 2 samples') + xlab('Samples from the same patient') + 
#   scale_fill_brewer(breaks = c('R3_R4','TI1_TI2','R4_TI2'), 
#                     labels = c('R inflamed x non-inflamed','TI inflamed x non-inflamed','R x TI non-inflamed '), 
#                     name = 'Locations compared')
# x_sig = wilcox.test(x = cor_df4$r[cor_df4$same_pnID  == F], y=cor_df4$r[cor_df4$same_pnID  == T])
# g4 = add_manual_significant_asterix( g4, val = cor_df4$r, lbl = x_sig$p.value, start_pos = 1, end_pos = 2, scale_var = 1.2, 
#                                      astx_pos = -1,asterix_scale_var = 1.2, label_size = 6, label_num_to_ast_flag = T )
# g1_sig = test_multiple_BH(var = cor_df4$xy_location_inflamation[cor_df4$same_pnID  == F], val = cor_df4$r[cor_df4$same_pnID  == F])
# 
# wilcox.test(x=cor_df4$r[cor_df4$same_pnID  == F & cor_df4$xy_location_inflamation == 'TI1_TI2'], 
#             y=cor_df4$r[cor_df4$same_pnID  == F & cor_df4$xy_location_inflamation == 'R3_R4'])
# out_fig = sprintf('%s/2_patients_taxa_bdiv_withtinBiopsy_chosen_location_inflamation_grouped%s.tiff',out_path, add_name)
# ggsave(out_fig,plot = g4, device = 'tiff', width = 6,height = 4, compression  = 'lzw', dpi = 300)
# 
# 
