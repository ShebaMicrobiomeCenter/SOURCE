library(ggplot2)
library(ggside)

source('/pita/users/tzipi/code/R_figs/figs_funcs.R')

side_type = 'density'
# side_type = 'boxplot'

# met_cols = c('SampleID','Dx','Patient_group2')
met_cols = c('SampleID','Dx','Patient_group2','location')

# path = '../res_china_stool_rural_ageMatched_amnonFlt/'
# name2 = 'rural_ageMatched'; cols = c('#7FB77E','#F7A76C')
# path = '../res_china_stool_rural_amnonFlt/'
# name2 = 'rural'; cols = c('#7FB77E','#F7A76C')

# path = '../res_israel_noRS_amnonFlt/'
# name2 = 'Israel'; cols = c('red','blue')
path = '../res_china_stool_urban_and_CD_amnonFlt/'
name2 = 'China_urban'; cols = c('red','blue')


name = 'source'
# taxa_file = sprintf('%s/biom/%s_deblur_map_L7.1.txt',path, name)
taxa_file = sprintf('%s/map.txt',path)
out_path =  sprintf('%s/R_res',path)
pcoa_file = sprintf('%s/core-metrics-results/unweighted_unifrac_pcoa_results/ordination.txt',path)
type = 'unweighted_unifrac'
na_str = c('no_data','_','NA','unknown', 'other','na','No_followup')
taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str)


out_path = sprintf('%s/PCoA',out_path)
dir.create(out_path)

res = pcoa_set_data(pcoa_file, taxa_file, met_cols)
pcoa = res[[1]]
dims = res[[2]]

# pcoa$V2 = -1*pcoa$V2
pcoa$Patient_group2[pcoa$Patient_group2=='Rural_health_<50%_in_city'] = 'Rural'
pcoa$Patient_group2[pcoa$Patient_group2=='Rural_health_>50%_in_city'] = 'Rural-Urban'
pcoa$Patient_group2[pcoa$Patient_group2=='Urban_health'] = 'Urban control'
pcoa$Patient_group2[pcoa$Patient_group2=='Chinese_crohns'] = 'CD'
pcoa$Patient_group2[pcoa$Patient_group2=='Israeli_Crohns'] = 'CD'
pcoa$Patient_group2[pcoa$Patient_group2=='Israeli_healthy'] = 'Control'

dim_labs = dims
for ( i in 1:3 )
  dim_labs[[i]] = sprintf("PC%d (%.2f%%)", i, dims[i]*100)

col_parm = 'Patient_group2'

pcoa$loc = ifelse(grepl('Stool',pcoa$location),'Stool','Biopsy') 

library(ggside)
g = ggplot(pcoa, aes_string(x='V2', y='V3' ,colour = col_parm)) + 
  geom_point(size = 2) + 
  # geom_point(size = 2, aes(shape = loc)) + scale_shape_manual(name = '', values = c(16,17)) + 
  scale_colour_manual(values = cols, name = '') +
  scale_fill_manual(values = cols, name = '') +
  theme_bw() + xlab(dim_labs[[1]]) + ylab(dim_labs[[2]]) +
  guides(alpha = "none") +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        text = element_text(size=16))
if (side_type == 'density')
{
  g = g+ geom_xsidedensity(aes(fill = Patient_group2, col = Patient_group2, alpha = 0.5)) +
    geom_ysidedensity(aes(fill = Patient_group2, col = Patient_group2, alpha = 0.5)) +
    theme(ggside.panel.scale = 0.25)
} else if (side_type == 'boxplot')
{
  g = g + geom_xsideboxplot(orientation = "y", aes(col = Patient_group2)) +
    geom_ysideboxplot(orientation = "x", aes(col = Patient_group2)) +
    theme(ggside.panel.scale = 0.25)
}

g
out_path = 'art_plots/res/'
# # ggsave(sprintf('%s/%s_pcoa_density.tiff', out_path, name2),
# # ggsave(sprintf('%s/%s_pcoa_boxplot.tiff', out_path, name2),
ggsave(sprintf('%s/%s_pcoa_%s.tiff', out_path, name2, side_type),
#        # plot = g, device = 'tiff', width =5.5,height = 3, compression='lzw')
       plot = g, device = 'tiff', width =5,height = 3, compression='lzw')

pc1_fig = ggplot(pcoa) +
  geom_boxplot(aes_string(x=col_parm, y='V2' ,fill = col_parm)) + # labs(colour = 'Day care') +
  theme_bw()  + ylab(dim_labs[[1]]) + 
  scale_fill_manual(values = cols, name = '') +
  theme(panel.grid = element_blank(), text = element_text(size=16), 
        axis.text.x = element_text(angle=45, hjust = 1))
res = add_significant_asterix_to_plot_BH_v2(p = pc1_fig, var = as.factor(pcoa[[col_parm]]), val = pcoa$V2, test_type = 'wilcox')
pc1_fig = res[[1]] 
# ggsave(sprintf('%s/%s_%s_PCoA_PC1_boxplot_wilcox_p%s.tiff', out_path, name2, col_parm, res[[2]]$p_val),plot = pc1_fig, device = 'tiff', width = 4,height = 3, compression  = 'lzw')
# # 
pc2_fig = ggplot(pcoa) +
  geom_boxplot(aes_string(x=col_parm, y='V3' ,fill = col_parm)) + # labs(colour = 'Day care') +
  theme_bw()  + ylab(dim_labs[[2]]) + 
  scale_fill_manual(values = cols, name = '') +
  theme(panel.grid = element_blank(), text = element_text(size=16), 
        axis.text.x = element_text(angle=45, hjust = 1))
res = add_significant_asterix_to_plot_BH_v2(p = pc2_fig, var = as.factor(pcoa[[col_parm]]), val = pcoa$V3, test_type = 'wilcox')
pc2_fig = res[[1]]
# ggsave(sprintf('%s/%s_%s_PCoA_PC2_boxplot_wilcox_p%s.tiff', out_path, name2, col_parm, res[[2]]$p_val),plot = pc2_fig, device = 'tiff', width = 4,height = 3, compression  = 'lzw')

