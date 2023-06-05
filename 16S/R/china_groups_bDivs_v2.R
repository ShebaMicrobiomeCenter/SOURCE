source('/pita/users/tzipi/code/R_figs/maaslin2_funs.R')
source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
library(Maaslin2)
library(ggplot2)
library(patchwork)

taxonomy_file = '/pita/users/tzipi/projects/16S/all_merges/DB1-20_merge/res/16S_DB1-20_merged/taxonomy/metadata_v2_source.tsv'

name = 'source'
path = '../res_china_stool_amnonFlt/'
# path = '../res_israel_noRS_amnonFlt/'
taxa_file = sprintf('%s/biom/%s_deblur_map_SVs.1.txt',path, name)
bdiv_file = sprintf('%s/core-metrics-results/unweighted_unifrac_distance_matrix/distance-matrix.tsv',path)
na_str = c('no_data','_','NA','unknown', 'other','na')
taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str)

bdiv = read.table(bdiv_file, header = T, row.names = 1, sep = '\t')

sv_df = taxa[,(which(names(taxa)=='Description')+1):dim(taxa)[2] ]
row.names(sv_df) = taxa$sample_ID
sv_df = as.data.frame(t(sv_df))

# sv_df = taxa[,c(1,(which(names(taxa)=='Description')+1):dim(taxa)[2]) ]
# sv_df = reshape2::melt(sv_df, value.name = 'RA')
# names(sv_df)[2] = 'SV'
# 
# 
# sv_df = merge(sv_df, taxa[,c('SampleID','pn_ID','Cohort','Patient_group','Dx','Location_inflamation','location')])

# dist_type = 'bDiv'; dist_lab = 'Unweighted Unifrac distance\nbetween 2 samples'
dist_type = 'corr'; dist_lab = 'Spearman correlations r\nof all taxa between 2 samples'

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

temp = data.frame(row.names = taxa$sample_ID, group = taxa$Patient_group2, SampleID = taxa$SampleID)
temp$group[temp$group == 'Chinese_crohns'] = 'CD'
temp$group[temp$group == 'Urban_health'] = 'Urban'
temp$group[temp$group == 'Rural_health_>50%_in_city'] = 'Rural-Urban'
temp$group[temp$group == 'Rural_health_<50%_in_city'] = 'Rural'
temp_x = temp[cor_df$x_pnID,]
temp_y = temp[cor_df$y_pnID,]

cor_df$x_group = temp_x$group
cor_df$y_group = temp_y$group
cor_df$Same_group = ifelse(cor_df$x_group==cor_df$y_group, 'Same group','Different groups')
cor_df$groups = sprintf('%s vs %s', cor_df$x_group, cor_df$y_group)

cor_df2 = cor_df
cor_df2 = cor_df2[cor_df2$Same_group == 'Different groups',]
cor_df2 = cor_df2[grepl('CD',cor_df2$groups) ,]

g = ggplot(cor_df2) + 
  geom_boxplot(aes(y=r, x=groups, fill = groups)) + 
  theme_classic() + ylab(dist_lab) + 
  # xlab('Samples from the same patient')+ 
  xlab('') + theme(legend.position = 'none') + 
  scale_fill_manual(values = c('gray80','gray80','gray80'), name = '')
# res_g = boxplot_2_vars_grouped_add_asterix( g, cor_df2, 
#                                             x_Val ='Same_patient', y_val = 'r', fill_val ='locations2', 
#                                             test_type = 'wilcox', fdr= 3)
test_type = 'perm'
# test_type = 'wilcox'
res_g = add_significant_asterix_to_plot_BH_v2(p = g, var = as.factor(cor_df2$groups), 
                                              val = cor_df2$r, 
                                              print_pvals = T,
                                              test_type = test_type)
out_path = 'res/group_correlations/'
# dir.create(out_path)
out_fig = sprintf('%s/%s_china_group_dist_boxplot_%s.tiff',out_path, dist_type, test_type)
res_g = format_fig(res_g[[1]], bw_flag = F, set_ttl_flag = F, x_angle = 45)
ggsave(out_fig,plot = res_g, device = 'tiff', width = 5,height = 5, 
       compression  = 'lzw', dpi = 300)


# g2 = ggplot(cor_df2, aes(x=locations, y=r)) + 
#   geom_boxplot(fill =  'gray') + 
#   facet_wrap(~Same_patient, scales = 'free_x') + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), 
#         panel.grid = element_blank()) + 
#   ggsignif::geom_signif( comparisons = list(c("Stool x R", "Stool x TI"),
#                                             c("Stool x R", "TI x R"),
#                                             c("Stool x TI", "TI x R")),
#                          step_increase = 0.1, map_signif_level = T, test = 'wilcox.test') + 
#   coord_cartesian(ylim = c(0.2,1)) + xlab('') + 
#   ylab(dist_lab) 


