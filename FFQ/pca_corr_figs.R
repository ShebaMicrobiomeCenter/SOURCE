source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
library(ggplot2)


col_parms = c('Patient_group2','Gender','Energy_Kcal_day')
col_parm = 'Patient_group2'

ffq_file = '../metadata/FFQ/FFQ_v11.txt'
ffq = read.table(ffq_file, header = T, sep = '\t', row.names = 1)

ffq = ffq[, !grepl('__n', names(ffq))]

food_prms = names(ffq)

map_file = '/pita/users/tzipi/projects/multiomics/SOURCE/metadata/alona_data/SOURCE_Israel_China_data_v11.txt'
map = read.table(map_file, header = T, sep = '\t')

metadata_prms = c('pn_ID','Patient_group','Time_in_City_last_year','Patient_group2','Age_years','BMI','Gender')

# map_f = merge(x = data.frame(pn_ID = row.names(ffq)), y = map[,metadata_prms], by.x = T )
map_f = map[map$pn_ID %in% row.names(ffq),]


# ## B300 looks bad! 41kg, 181cm, 7k calories per day? outlier to remove
# df = df[df$pn_ID != 'B300',]

# name = 'China'
# groups = c('Rural_health_<50%_in_city','Rural_health_>50%_in_city','Urban_health','Chinese_crohns')
# cols = c('#82CD47','#1F8A70','#3DB2FF','#EB5353')
# groups2 = c('Rural','Rural-Urban','Urban','CD')

# name = 'China_rural'
# groups = c('Rural_health_<50%_in_city','Rural_health_>50%_in_city')
# cols = c('#82CD47','#1F8A70','#3DB2FF','#EB5353')

name = 'China_urban'
groups = c('Urban_health','Chinese_crohns')
cols = c('#3DB2FF','#EB5353')


# name = 'Israel'
# groups = c('Israeli_healthy','Israeli_Crohns')
# cols = c('#3DB2FF','#EB5353')
# groups2 = c('Control','CD')


cor_type = 'spearman'
# cor_type = 'pearson'

map_f = map_f[map_f$Patient_group2 %in% c(groups) & !is.na(map_f$Patient_group2),]
ffq = ffq[map_f$pn_ID,]


out_path = 'res_v11/pca_biplot/'
dir.create(out_path)

df2 = ffq
# row.names(df2) = df$pn_ID

good_pos = c()
for ( i in 1:dim(df2)[2] )
{
  if ( sum(!is.na(df2[,i]) ) > 5 & sum(!df2[,i] ==0, na.rm = T )!=0 )
    # if ( sum(is.na(df2[,i]) ) == 0  & sum(!df2[,i] ==0, na.rm = T )!=0 )
  {
    good_pos = c(good_pos, i)
  }
}
df2 = df2[,good_pos]

good_pos = c()
for ( i in 1:dim(df2)[1] )
{
  if ( sum(!is.na(df2[i,]), na.rm = T) != 0  )
    # if ( sum(is.na(df2[i,])) == 0  )
  {
    good_pos = c(good_pos, i)
  }
}
df2 = df2[good_pos,]

df2 = df2[ rowSums( is.na(df2) ) == 0, ]

map_f = map_f[map_f$pn_ID %in% row.names(df2),]

df3 = df2
# for ( i in 1:dim(df3)[2] )
# {
#   df3[,i] = (df3[,i] - mean( df3[,i], na.rm=T ) ) / sd(df3[,i], na.rm = T)
#   # df3[,i] = (df3[,i] / mean(df3[,i]))
# }

df.pca = prcomp(df3, center = T, scale = T)
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
# ggsave(sprintf('%s/%s_%s_%s_pca.tiff', out_path, paste(groups, collapse = '_'), col_parm, type ),plot = g, device = 'tiff', width = 5,height =2.5, compression  = 'lzw')
# ggsave(sprintf('%s/%s_%s_pca.tiff', out_path, paste(groups, collapse = '_'), col_parm ),
#        plot = g, device = 'tiff', width =5.5,height = 2.5, compression  = 'lzw')
ggsave(sprintf('%s/%s_%s_pca_density.tiff', out_path, paste(groups, collapse = '_'), col_parm ),
# ggsave(sprintf('%s/%s_%s_pca_boxplot.tiff', out_path, paste(groups, collapse = '_'), col_parm ),
       plot = g, device = 'tiff', width =5.5,height = 3, compression  = 'lzw')
       # plot = g, device = 'tiff', width =5,height = 3, compression  = 'lzw')

temp = ggplot(data) + geom_boxplot(aes(x=map_f$Patient_group2, y=data$PC1))
res = add_significant_asterix_to_plot_BH_v2(p = temp, 
                                            var = as.factor(map_f$Patient_group2),
                                            val = data$PC1)
pc1_fig = ggplot(data) +
  geom_boxplot(aes_string(x=map_f$Patient_group2, y=data$PC1 ,fill = map_f$Patient_group2)) + # labs(colour = 'Day care') +
  theme_bw()  + ylab(dim_labs[[1]]) + xlab('') + 
  scale_fill_manual(values = cols, name = '') +
  theme(panel.grid = element_blank(), text = element_text(size=16), 
        axis.text.x = element_text(angle=45, hjust = 1), legend.position = 'none')
res = add_significant_asterix_to_plot_BH_v2(p = pc1_fig, var = as.factor(map_f$Patient_group2), 
                                            val = data$PC1, test_type = 'wilcox', 
                                            show_pval_as_asterix = F, 
                                            label_size = 3, asterix_scale_var = 2.3)
pc1_fig = res[[1]] 
ggsave(sprintf('%s/%s_%s_pc1_boxplot.tiff', out_path, paste(groups, collapse = '_'), col_parm ),
        plot = pc1_fig, device = 'tiff', width =4,height = 4, compression  = 'lzw')

pc2_fig = ggplot(data) +
  geom_boxplot(aes_string(x=map_f$Patient_group2, y=data$PC2 ,fill = map_f$Patient_group2)) + # labs(colour = 'Day care') +
  theme_bw()  + ylab(dim_labs[[2]]) + xlab('') + 
  scale_fill_manual(values = cols, name = '') +
  theme(panel.grid = element_blank(), text = element_text(size=16), 
        axis.text.x = element_text(angle=45, hjust = 1), legend.position = 'none')
res = add_significant_asterix_to_plot_BH_v2(p = pc2_fig, var = as.factor(map_f$Patient_group2), 
                                            val = data$PC2, test_type = 'wilcox', 
                                            show_pval_as_asterix = F, 
                                            label_size = 3, asterix_scale_var = 2.3)
pc1_fig = res[[1]] 
ggsave(sprintf('%s/%s_%s_pc2_boxplot.tiff', out_path, paste(groups, collapse = '_'), col_parm ),
       plot = pc1_fig, device = 'tiff', width =4,height = 4, compression  = 'lzw')




mets = df3
cor_r_p1 = vector(mode = 'numeric',length = dim(mets)[2] )
cor_p_p1 = vector(mode = 'numeric',length = dim(mets)[2] )
cor_r_p2 = vector(mode = 'numeric',length = dim(mets)[2] )
cor_p_p2 = vector(mode = 'numeric',length = dim(mets)[2] )
for ( i in 1:dim(mets)[2] )
{
  res_p1 = cor.test(x = data$PC1, y=mets[,i], method = 'spearman')
  res_p2 = cor.test(x = data$PC2, y=mets[,i], method = 'spearman')
  cor_r_p1[i] = res_p1$estimate
  cor_r_p2[i] = res_p2$estimate
  cor_p_p1[i] = res_p1$p.value
  cor_p_p2[i] = res_p2$p.value
}
df = data.frame(cor_r_p1 = cor_r_p1, cor_r_p2 = cor_r_p2, cor_p_p1 = cor_p_p1, cor_p_p2 = cor_p_p2, parameter = names(mets))
# df = df[ ( df$cor_p_p1 <= 0.05 & abs(df$cor_r_p1) > 0.5)  |  ( df$cor_p_p2 <= 0.05 & abs(df$cor_r_p2) > 0.5), ]
df = df[ df$cor_p_p1 <= 0.05 | df$cor_p_p2 <= 0.05, ]


unique_cols = c('#e6194b', '#3cb44b', '#ffe119',
                '#4363d8', '#f58231', '#911eb4', '#46f0f0',
                '#f032e6', '#bcf60c', '#9a6324', '#008080',
                '#e6beff', '#aaffc3', '#800000', '#fabebe',
                '#fffac8', '#808000', '#ffd8b1', '#000075',
                '#808080', '#000000') # '#ffffff' = white
df$r_sums = max(abs(df$cor_r_p1), abs(df$cor_r_p2))
# df$r_sums = abs(df$cor_r_p1) + abs(df$cor_r_p2)
df = df[order(df$r_sums, decreasing = T),]
df$parameter = factor(df$parameter, levels = df$parameter)

df = df[1:10,]

gc2 = ggplot(df) + 
  geom_segment(data = df, aes(x = 0, y = 0, xend = (cor_r_p1), yend = (cor_r_p2)),
               size = 0.5, color = scales::muted('blue'),alpha = 0.7,
               arrow = arrow(length = unit(0.5, "picas"))) + 
  ggrepel::geom_text_repel(data=df, aes(x=cor_r_p1, y=cor_r_p2, label=parameter), 
                  size = 2, color=scales::muted('blue')) +
  theme_bw() + 
  xlab('PC1 corelation R') + ylab('PC2 corelation R') + 
  scale_colour_manual(name = '',values = unique_cols) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #  + 
ggsave(sprintf('%s/%s_%s_pca_corrs2.tiff', out_path, paste(groups, collapse = '_'), col_parm ),
       plot = gc2, device = 'tiff', width =3.5,height = 2.5, compression  = 'lzw')
# geom_point(aes(x=cor_r_p1, y=cor_r_p2,colour = parameter)) # + 
# geom_text(data = df, aes(label = parameter, x = (cor_r_p1*(4/5)), y = (cor_r_p2)*(4/5)), 
# geom_text(data = df, aes(label = parameter, x = (cor_r_p1), y = (cor_r_p2)), 
#           color = "black", size = 4,  angle = df$Angle)#, vjust = df$Offset )


# ggsave(sprintf('%s/PCA_corrs_%s_%s.pdf', out_path,c,type),plot = g, device = 'pdf', width = 5,height = 3)


# gggsave(sprintf('%s/%s_%s_biplot.tiff', out_path, paste(groups, collapse = '_'), col_parm ),plot = gb, device = 'tiff', width =8,height = 5, compression  = 'lzw')

# }
