
library(ggplot2)
library(stringr)

source('/pita/users/tzipi/code/R_figs/figs_funcs.R')

source('/pita/users/tzipi/code/R_figs/maaslin_funcs.R')

metadata_cols = c('SampleID','dysbiosis_index','faith_pd',
                  'Patient_group2','pn_ID')

check_val = 'Patient_group2'
na_str = c('no_data','_','NA','unknown', 'other','na','No_followup','ND')

name = 'source'
if ( name == 'source')
{
  path = '../res_china_stool_rural_ageMatched_amnonFlt/'
  # path = '../res_china_stool_rural_amnonFlt/'
  taxa_file = sprintf('%s/biom/%s_deblur_map_L7.1.txt',path, name)
  SV_file = sprintf('%s/biom/%s_deblur_map_SVs.1.txt',path, name)
  a_file = sprintf('%s/core-metrics-results/faith_pd_vector/alpha-diversity.tsv',path)
  b_file = sprintf('%s/core-metrics-results/unweighted_unifrac_distance_matrix/distance-matrix.tsv',path)
  pcoa_file = sprintf('%s/core-metrics-results/unweighted_unifrac_pcoa_results/ordination.txt',path)
  out_path =  sprintf('%s/R_res',path)
}
dir.create(out_path)

# taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str)
# map = taxa[,1:which(names(taxa)=='Description') ]
# map = map[map$Time==1,]
# write.table(x = map,file = sprintf('%s/DB_2_3_5_7_10_11_14_lungs_ctrl_map_v8_2_time1.txt','../data/'),quote = F,sep = '\t',row.names = F)

bdiv_dist_file = sprintf('%s/core-metrics-results/unweighted-unifrac-Patient_group2-significance/raw_data.tsv',path)
bdiv = read.table(bdiv_dist_file,sep="\t", header=TRUE, na.strings = na_str)

bdiv$Dx_res = sprintf('%s-%s',bdiv$Group1, bdiv$Group2)
bdiv = bdiv[bdiv$Group1 == bdiv$Group2,]

bdiv$Dx_res = gsub('.*-','',bdiv$Dx_res )
bdiv$Dx_res = gsub('_',' ',bdiv$Dx_res )
# bdiv$Dx_res = gsub(' health','',bdiv$Dx_res )
# bdiv$Dx_res = gsub(' in city',' time\nin city',bdiv$Dx_res )

bdiv$Dx_res[bdiv$Dx_res == 'Rural health <50% in city'] = 'Rural'
bdiv$Dx_res[bdiv$Dx_res == 'Rural health >50% in city'] = 'Rural-Urban'

g = ggplot(bdiv) + geom_boxplot(aes(x=Dx_res, y=Distance), fill = 'gray75') + 
  # scale_fill_manual(values=c('#0099FF','purple','#FF6666')) + 
  theme_bw() + ylab('Unweighted Unifrac Distance' ) + xlab(check_val) 
g = format_fig(g,set_ttl_flag = F, x_angle = 45)+ theme(legend.position = "none")
# ggsave(sprintf('%s/%s_time1.tiff', out_path, 'beta_boxplot'),plot = g, device = 'tiff', width =5,height = 5, compression='lzw')
# res = add_significant_asterix_to_plot_BH(p = g, var = as.factor(bdiv$Dx_res), val = bdiv$Distance, print_pvals = T)

# ggsave(sprintf('%s/%s_time1.tiff', out_path, 'beta_boxplot'),plot = res[[1]], device = 'tiff', width =5,height = 5, compression='lzw')

bdiv_sig_file = sprintf('%s/core-metrics-results/unweighted-unifrac-Patient_group2-significance/permanova-pairwise.csv',path)
bdiv_sig = read.csv2(bdiv_sig_file,sep=',', header=TRUE, na.strings = na_str)

out_path = 'art_plots/res/'
g = add_manual_significant_asterix(g, val = bdiv$Distance, lbl = as.numeric(bdiv_sig$q.value), 
                                   label_size = 8, label_num_to_ast_flag = T, scale_var = 0.5 ) + 
  xlab('') + theme(panel.grid = element_blank())
ggsave(sprintf('%s/%s_%s_qiime2_perm_test.tiff', out_path, 'beta_boxplot', check_val),
#        plot = g, device = 'pdf', width =2,height = 3.5)
       plot = g, device = 'tiff', width =2,height = 3.5, compression = 'lzw',dpi=600)


