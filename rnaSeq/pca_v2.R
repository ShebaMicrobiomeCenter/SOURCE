source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
library(ggplot2)
library(ggside)

out_path = 'pca_v2/'
dir.create(out_path)

# country = 'Israel'
country = 'China'

## read and organize data
name = 'SOURCE_israel_china'
tpm = read.table(file = sprintf('data/%s_rnaSeq_txi_res.txt',name),header = T,sep = '\t',row.names = 1)
map = read.table('data/SOURCE_israel_china_rnaSeq_map.txt',header = T,sep = '\t')
map$SampleID = make.names(map$SampleID)


df = tpm


col_parms = c('Patient_group2')
col_parm = 'Patient_group2'


name = 'China'
groups = c('Rural_health_<50%_in_city','Rural_health_>50%_in_city','Urban_health','Chinese_crohns')
# cols = c('#82CD47','#1F8A70','#3DB2FF','#EB5353')
cols = c('#82CD47','#3DB2FF','#EB5353')
groups2 = c('Rural','Rural-Urban','Urban','CD')

# name = 'China_rural'
# groups = c('Rural_health_<50%_in_city','Rural_health_>50%_in_city')
# cols = c('#82CD47','#1F8A70','#3DB2FF','#EB5353')

# name = 'China_urban'
# groups = c('Urban_health','Chinese_crohns')
# cols = c('#3DB2FF','#EB5353')


# name = 'Israel'
# groups = c('Israeli_healthy','Israeli_Crohns')
# cols = c('#3DB2FF','#EB5353')
# groups2 = c('Control','CD')


cor_type = 'spearman'
# cor_type = 'pearson'

map_f = map
map_f = map_f[map_f$Patient_group2 %in% c(groups) & !is.na(map_f$Patient_group2),]
tpm = tpm[,map_f$SampleID]


# out_path = 'res_v11/pca_biplot/'
# dir.create(out_path)

df = tpm
# row.names(df2) = df$pn_ID
# filtering to wanted genes (protein coding from the table yael send
gene_file = '/pita/users/tzipi/projects/rnaSeq/lncRNA/tables/from_yael/19815_Protein_coding_genes_GencodeV24.txt'
gene_df = read.table(gene_file, header = T)
wanted_genes = sprintf('%s_%s',gene_df$Gene, gene_df$Ensembl_Gene_ID)
txi_genes = row.names(df)
wanted_genes_pos1 = txi_genes %in% wanted_genes

## filter TPM to 1 TPM>0.2%
tpm_genes = row.names(df)
wanted_genes_pos2 = vector(mode = 'logical',length = length(tpm_genes))
for ( i in 1:length(tpm_genes) )
{
  per = sum(df[i,] > 1) / dim(df)[2]
  wanted_genes_pos2[i] = per > 0.2
}
wanted_genes_pos = wanted_genes_pos1 & wanted_genes_pos2

df = df[wanted_genes_pos, ]

## organize data for pca analysis
tpm = df
tpm2 = tpm[rowSums(is.na(tpm)) == 0,] # remove any line with >1 NA
tpm2 = tpm2[rowSums(tpm2) != 0,] # remove any gene with 0 expression
input_data = tpm2


df3 = tpm
df.pca = prcomp(t(df3), center = T, scale = T)
imp = summary(df.pca)$importance[2,]
data = as.data.frame(df.pca$x)
# col_parm = 'Patient_group'

# for(col_parm in col_parms)
# {

map_f$Patient_group2[map_f$Patient_group2=='Rural_health_<50%_in_city'] = 'Rural'
map_f$Patient_group2[map_f$Patient_group2=='Rural_health_>50%_in_city'] = 'Rural-Urban'
map_f$Patient_group2[map_f$Patient_group2=='Urban_health'] = 'Urban'
map_f$Patient_group2[map_f$Patient_group2=='Chinese_crohns'] = 'CD'
map_f$Patient_group2[map_f$Patient_group2=='Israeli_Crohns'] = 'CD'
map_f$Patient_group2[map_f$Patient_group2=='Israeli_healthy'] = 'Control'


map_f$Patient_group2 = factor(map_f$Patient_group2, levels = groups2)
data$Patient_group2 = map_f$Patient_group2
# df$Patient_group2 = factor(df$Patient_group2, levels = c('Rural_health_<50%_in_city','Rural_health_>50%_in_city','Urban_health','Chinese_crohns'))
g = ggplot(data, aes(x=PC1, y=PC2)) +
  geom_point(aes(colour = map_f[[col_parm]] ), size=2) + 
  # scale_colour_manual(values =c('#3cb44b','#ffe119','#f58231','#9b2226')) +
  xlab(sprintf('PC1 (%.2f%%)', imp[1]*100)) + ylab(sprintf('PC2 (%.2f%%)', imp[2]*100)) +
  # guides(colour=guide_legend(title=NULL)) +
  scale_color_manual(values = cols, name='') +
  scale_fill_manual(values = cols, name='') +
  geom_xsidedensity(aes(fill = data[[col_parm]], col = data[[col_parm]], alpha = 0.5)) +
  geom_ysidedensity(aes(fill = data[[col_parm]], col = data[[col_parm]], alpha = 0.5)) +
  # geom_xsideboxplot(orientation = "y", aes(col = data[[col_parm]]), position = 'dodge2') +
  # geom_ysideboxplot(orientation = "x", aes(col = data[[col_parm]])) +
  theme_bw() + 
  guides(alpha = "none") + 
  theme(panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        text = element_text(size=16), ggside.panel.scale = 0.25)
g
ggsave(sprintf('%s/%s_%s_pca_density.tiff', out_path, paste(groups, collapse = '_'), col_parm ),
       # ggsave(sprintf('%s/%s_%s_pca_boxplot.tiff', out_path, paste(groups, collapse = '_'), col_parm ),
       plot = g, device = 'tiff', width =5.5,height = 3, compression  = 'lzw')
# plot = g, device = 'tiff', width =5,height = 3, compression  = 'lzw')

map_f = droplevels(map_f)
pc1_fig = ggplot(data) +
  geom_boxplot(aes_string(x=map_f$Patient_group2, y=data$PC1 ,fill = map_f$Patient_group2)) + # labs(colour = 'Day care') +
  theme_bw()  +  ylab(sprintf('PC1 (%.2f%%)', imp[1]*100)) + xlab('') + 
  scale_fill_manual(values = cols, name = '') +
  theme(panel.grid = element_blank(), text = element_text(size=16), 
        axis.text.x = element_text(angle=45, hjust = 1), legend.position = 'none')
res = add_significant_asterix_to_plot_BH_v2(p = pc1_fig, var = as.factor(map_f$Patient_group2), 
                                            val = data$PC1, test_type = 'wilcox', 
                                            show_pval_as_asterix = T, print_pvals = T
                                            # label_size = 3, asterix_scale_var = 2.3
                                            )
pc1_fig = res[[1]] 
ggsave(sprintf('%s/%s_%s_pc1_boxplot.tiff', out_path, paste(groups, collapse = '_'), col_parm ),
       plot = pc1_fig, device = 'tiff', width =4,height = 4, compression  = 'lzw')

pc2_fig = ggplot(data) +
  geom_boxplot(aes_string(x=map_f$Patient_group2, y=data$PC2 ,fill = map_f$Patient_group2)) + # labs(colour = 'Day care') +
  theme_bw()  + ylab(sprintf('PC2 (%.2f%%)', imp[2]*100)) + xlab('') + 
  scale_fill_manual(values = cols, name = '') +
  theme(panel.grid = element_blank(), text = element_text(size=16), 
        axis.text.x = element_text(angle=45, hjust = 1), legend.position = 'none')
res = add_significant_asterix_to_plot_BH_v2(p = pc2_fig, var = as.factor(map_f$Patient_group2), 
                                            val = data$PC2, test_type = 'wilcox', 
                                            show_pval_as_asterix = T, print_pvals = T
                                            # label_size = 3, asterix_scale_var = 2.3
                                            )
pc2_fig = res[[1]] 
ggsave(sprintf('%s/%s_%s_pc2_boxplot.tiff', out_path, paste(groups, collapse = '_'), col_parm ),
       plot = pc2_fig, device = 'tiff', width =4,height = 4, compression  = 'lzw')

