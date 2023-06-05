source('/pita/users/tzipi/code/R_figs/figs_funcs.R')

type = 'v11'
ffq_file = '../metadata/FFQ/FFQ_v11.txt'
ffq = read.table(ffq_file, header = T, sep = '\t', row.names = 1)
ffq = ffq[, !grepl('__n', names(ffq))]

# clust_type = 'ward.d'
clust_type = 'heatmapClust'

food_prms = names(ffq)

map_file = '/pita/users/tzipi/projects/multiomics/SOURCE/metadata/alona_data/SOURCE_Israel_China_data_v11.txt'
map = read.table(map_file, header = T, sep = '\t')

metadata_prms = c('pn_ID','Patient_group','Time_in_City_last_year','Patient_group2','Age_years','BMI','Gender')

# map_f = merge(x = data.frame(pn_ID = row.names(ffq)), y = map[,metadata_prms], by.x = T )
map_f = map[map$pn_ID %in% row.names(ffq),]


# ## B300 looks bad! 41kg, 181cm, 7k calories per day? outlier to remove
# df = df[df$pn_ID != 'B300',]

# groups = c('Rural_health')
# groups = c('Urban_health','Chinese_crohns')
groups = c('Urban_health','Rural_health','Chinese_crohns')
# groups = c('Israeli_healthy')
# groups = c('Israeli_healthy','Israeli_Crohns')
# cor_type = 'spearman'
# cor_type = 'pearson'

map_f = map_f[map_f$Patient_group %in% c(groups) & !is.na(map_f$Patient_group),]
ffq = ffq[map_f$pn_ID,]

out_path = 'res_v11/corr_heatamaps/'
dir.create(out_path)

df2 = ffq
# row.names(df2) = df$pn_ID

good_pos = c()
for ( i in 1:dim(df2)[2] )
{
  if ( sum(!is.na(df2[,i])) > 5 & sum(!df2[,i] ==0, na.rm = T )!=0 )
  {
    good_pos = c(good_pos, i)
  } else
    print(names(df2)[i])
}
df2 = df2[,good_pos]

good_pos = c()
for ( i in 1:dim(df2)[1] )
{
  if ( sum(!is.na(df2[i,])) > 2  )
  {
    good_pos = c(good_pos, i)
  }
}
df2 = df2[good_pos,]


# df3 = df2
# # df3 = log(df3+0.0001)
# for ( i in 1:dim(df3)[2] )
# {
#   df3[,i] = (df3[,i] - mean(df3[,i], na.rm = T)) / sd(df3[,i], na.rm = T) 
#   # df3[,i] = (df3[,i] / mean(df3[,i])) 
# }

# install.packages('PerformanceAnalytics')
# library("PerformanceAnalytics")
# chart.Correlation(df2, histogram=TRUE, pch=19)


# res <- cor(df2, method = cor_type, na.rm=T)
library("Hmisc")
res2 <- rcorr(as.matrix(df2), type = cor_type)


col2 = colorRampPalette(rev( c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                               "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061") ))(200)
# library(corrplot)
# corrplot(res2$r, type = "upper", order = "hclust", 
#          tl.col = "black", method = 'circle', col = col2, tl.cex= 0.8)

ht = heatmap(x = res2$r, col = col2, symm = T, scale = 'none')

library(reshape)
res_melt = reshape::melt(data = res2$r)
temp = reshape::melt(data = res2$P)
res_melt$p_val = temp$value

## add q values without repeats
temp = res_melt[!duplicated(t(apply(res_melt[,1:2], 1, sort))),]
temp$q_val = p.adjust(temp$p_val, method = 'BH', n = length(temp$p_val))
res_melt$q_val = -1
for ( i in 1:dim(res_melt[1]))
{
  pos = which( ( res_melt$X1[i] == temp$X1 & res_melt$X2[i] == temp$X2 ) |
                 ( res_melt$X1[i] == temp$X2 & res_melt$X2[i] == temp$X1 ) )
  if (length(pos)!=1)
    print('problem')
  res_melt$q_val[i] = temp$q_val[pos]
}

# data <- scale(t(res2$r))
# ord <- hclust( dist(data, method = "euclidean"), method = "ward.D" )$order
# ord <- hclust( dist(data, method = "euclidean"), method = "complete" )$order
if (clust_type == 'ward.d')
{
  ord = hclust( dist( res2$r, method = "euclidean"), method = "ward.D" )$order
  # ord = hclust( dist( res2$r))$order
} else if (clust_type == 'heatmapClust')
{
  ord = ht$rowInd
}

res_melt$X1 = factor(res_melt$X1, levels = row.names(res2$r)[ord])
res_melt$X2 = factor(res_melt$X2, levels = rev(row.names(res2$r)[ord]))

cor_heatmap = ggplot(res_melt, aes(x=X1, y=X2, fill = value)) + 
  geom_tile() + 
  # scale_fill_gradientn(colours = col2) +
  scale_fill_gradientn(colours = col2, limits = c(-1,1), name = 'Spearman\'s\nrho') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  # geom_point(aes(shape = q_val<=0.05), size=0.7) + scale_shape_manual(values = c(NA, 20)) +
  # geom_text(aes(label=sprintf('%.2f\n%.2e',value,q_val )), size=3)
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  # geom_text(aes(label=sprintf('%.2f',value )), size=3) + 
  xlab('') + ylab('')

# ggsave(sprintf('%s/cor_heatmap_%s_%s_%s_%s.pdf', out_path, type, cor_type, clust_type,
#                paste(groups, collapse = '_') ),plot = cor_heatmap,  device = 'pdf',
#                width = 7,height = 6)
# #                width = 10,height = 9)


# df3 = df2
# for ( i in 1:dim(df3)[2] )
# {
#   df3[,i] = (df3[,i] - mean(df3[,i])) / sd(df3[,i]) 
#   # df3[,i] = (df3[,i] / mean(df3[,i])) 
# }
# 
# 
# colpn2 = df$Patient_group2
# colpn2[colpn2 == 'Urban_health'] = '#f58231'
# colpn2[colpn2 == 'Rural_health_<50%_in_city'] = '#3cb44b'
# colpn2[colpn2 == 'Rural_health_>50%_in_city'] = '#ffe119'
# colpn2[colpn2 == 'Rural_health_10-50%_in_city'] = '#911eb4'
# 
# colGender = df$Gender
# colGender[colGender == 'male'] = '#42d4f4'
# colGender[colGender == 'female'] = '#f032e6'
# 
# Heatmap(t(as.matrix(df3)),  col = col2, ColSideColors  = colpn2 )
# 
# temp = df3; temp$SampleID = row.names(temp)
# z_melt = melt(data = temp, id.vars = 'SampleID')
# 
# ord <- hclust( dist(scale(t(df3)), method = "euclidean"), method = "ward.D" )$order
# z_melt$variable = factor(z_melt$variable, names(df3)[ord])
# ord <- hclust( dist(scale(df3), method = "euclidean"), method = "ward.D" )$order
# z_melt$SampleID = factor(z_melt$SampleID, row.names(df3)[ord])
# 
# z_melt2 = z_melt
# z_melt2$value[z_melt2$value > 5] = 5
# z_heatmap = ggplot(z_melt2, aes(x=SampleID, y=variable, fill = value)) + 
#   # geom_tile() + 
#   scale_fill_gradientn(colours = col2) +
#   # scale_fill_gradientn(colours = col2, limits = c(-1,1)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
#   # geom_point(aes(shape = q_val<=0.05)) + scale_shape_manual(values = c(NA, 8)) +
#   # geom_text(aes(label=sprintf('%.2f\n%.2e',value,q_val )), size=3)
#   # geom_text(aes(label=sprintf('%.2f',value )), size=3) + 
#   xlab('') + ylab('')
# # ggsave(sprintf('%s/z_score_heatmap_%s.tiff', out_path, paste(groups, collapse = '_') ),plot = z_heatmap, device = 'tiff', width = 25,height = 4, compression  = 'lzw')
# # # ggsave(sprintf('%s/z_score_heatmap_%s.tiff', out_path, paste(groups, collapse = '_') ),plot = z_heatmap, device = 'tiff', width = 8,height = 10, compression  = 'lzw')
# 
# unique_colors = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
# df.pca = prcomp(df3, center = T, scale = T)
# data = as.data.frame(df.pca$x)
# # col_parm = 'Patient_group'
# col_parm = 'Patient_group2'
# # col_parm = 'X16d._How_long_have_you_stayed_in_a_city_in_the_last_year_.'
# # col_parm = 'Fat_gr_day'
# g = ggplot(data, aes(x=PC1, y=PC2)) +
#   # geom_point(colour='gray') + 
#   # geom_point(aes(colour = df$Patient_group, shape = df$X16d._How_long_have_you_stayed_in_a_city_in_the_last_year_.)) + 
#   geom_point(aes(colour = df[[col_parm]] )) + guides(colour=guide_legend(title=col_parm)) +
#   scale_colour_manual(values =unique_colors) + 
#   # geom_point() + 
#   # geom_text(aes(label=row.names(data)), size=3)+
#   # geom_text(aes(label=names(tpm2), colour = map$loc), size=3) +
#   theme_bw() # + scale_colour_brewer(palette = 'Paired')
# g
# # ggsave(sprintf('%s/pca_%s_%s.tiff', out_path, col_parm, paste(groups, collapse = '_') ),plot = g, device = 'tiff', width = 6,height = 4, compression  = 'lzw')

## make fat iron scatter
if ( all(groups) %in% c('Urban_health','Rural_health','Chinese_crohns') )
{
  cols = c('#82CD47','#1F8A70','#3DB2FF','#EB5353')
  groups2 = c('Rural','Rural-Urban','Urban','CD')
  map_f$Patient_group2[map_f$Patient_group2=='Rural_health_<50%_in_city'] = 'Rural'
  map_f$Patient_group2[map_f$Patient_group2=='Rural_health_>50%_in_city'] = 'Rural-Urban'
  map_f$Patient_group2[map_f$Patient_group2=='Urban_health'] = 'Urban'
  map_f$Patient_group2[map_f$Patient_group2=='Chinese_crohns'] = 'CD'
  map_f$Patient_group2 = factor(map_f$Patient_group2, levels = groups2)
  g = ggplot(df2, aes(x=Ferrus_mg_day, y=Fat_gr_day, colour = map_f$Patient_group2)) + 
    geom_point() + theme_bw() + 
    scale_colour_manual(values = cols, name='') + 
    theme(panel.grid = element_blank())

  ggsave(sprintf('%s/../fat_iron_china_scatter.pdf', out_path),
         plot = g,  device = 'pdf',width = 5,height = 4)
}

