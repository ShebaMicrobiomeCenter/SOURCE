library(ggplot2)

get_good_enrichment_names_list = function(df, top_num = 30)
{
  # df = df[df$Category %in% c('GO: Molecular Function','GO: Biological Process','GO: Cellular Component',
  #                            'Pathway','Coexpression Atlas','ToppCell Atlas'),]
  good_names = c()
  cats = unique(df$Category)
  mds = unique(df$module)
  for (md in mds)
  {
    for( cat in cats )
    {
      temp = df[df$Category == cat & df$module == md & !is.na(df$FDR_BH),]
      temp = temp[base::order(temp$FDR_BH, decreasing =T),]
      temp = temp[1:top_num,]
      good_names = c(good_names,temp$Name)
    }
  }
  good_names = unique(good_names)
  return(good_names)
}
wanted_modules = c('yellow','green','red','pink','purple','tan','salmon','black','brown')
top_num = 10

# source('/pita/users/tzipi/projects/multiomics/SOURCE/multi_analysis/WGCNA/rnaSeq_main/modules_enrichment/org_function_data.R')
data_file = 'rnaSeq_main/modules_enrichment/data/SOURCE_Israel_TI_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default/enrichment_files/toppgene_merged_ens.txt'

df = read.table(data_file, sep = '\t',header = T, comment.char = '', quote="\"", na.strings = c('NA','na','',' '))
df_mat = df
# fdr_cut = 120
# good_pos = apply(df[,4:12], 1, function(x) max(x)) > fdr_cut
# df = df[good_pos,]

# df = df[df$Category %in% c('GO: Molecular Function','GO: Biological Process','GO: Cellular Component',
#                               'Pathway','Coexpression Atlas','ToppCell Atlas'),]
# df = df[df$Category %in% c('GO: Molecular Function','GO: Biological Process','GO: Cellular Component',
#                            'Pathway'),]

bad_pos = df$Category == 'Coexpression Atlas' & ! grepl(pattern = 'immune',x = tolower(df$Name))
df = df[-which(bad_pos),]
bad_pos = df$Category == 'ToppCell Atlas' & ! grepl(pattern = 'ileum',x = tolower(df$Name))
df = df[-which(bad_pos),]


df = reshape::melt(data = df, id.vars = c('Name','Category','ID'),variable_name = 'module')
names(df)[names(df) == 'value']  = 'FDR_BH'
df$FDR_BH[df$FDR_BH <=0.05] = NA


good_names = get_good_enrichment_names_list(df, top_num = top_num)
df = df[df$Name %in% good_names,]

# ## filter to top x functions
# top_num = 50
# temp = df[!is.na(df$FDR_BH),]
# temp = temp[base::order(temp$FDR_BH, decreasing =T),]
# df = df[df$Name %in% temp$Name[1:top_num],]

df = droplevels(df)

# if(all(temp$red==0))
# {
#   pos = which(names(temp)=='red')
#   temp = temp[,-pos]
#   df = df[df$module!='red',]
#   interesting_modules = c('blue','yellow','black','green', 'brown','pink')
# }

# heatmap(as.matrix(temp), scale = 'none')

temp = df_mat[df_mat$Name %in% df$Name,]
ord = hclust( dist(temp[,4:12], method = "euclidean") )$order
df$Name = factor( df$Name, levels = unique(df$Name)[ord] )

# ord = order( apply(temp[,c('blue','yellow')], 1, mean )  ) 
# dff$Name = factor( dff$Name, levels = sort(unique(dff$Name))[ord] )

# ord = order( apply(temp[,c('black','green','brown','pink')], 1, mean )  ) 
# dff$Name = factor( dff$Name, levels = sort(unique(dff$Name))[ord] )

# ord = order( apply(temp[,c('black','green','brown','pink')], 1, mean, na.rm=T ) - apply(temp[,c('blue','yellow')], 1, mean, na.rm=T ) ) 
# df$Name = factor( df$Name, levels = sort(unique(df$Name))[ord] )


# df$module = factor(df$module, levels = interesting_modules)
# 
# unique_cols = c('#e6194b', '#3cb44b', '#ffe119',
#                 '#4363d8', '#f58231', '#911eb4', '#46f0f0',
#                 '#f032e6', '#bcf60c', '#9a6324', '#008080',
#                 '#e6beff', '#aaffc3', '#800000', '#fabebe',
#                 '#fffac8', '#808000', '#ffd8b1', '#000075',
#                 '#808080', '#000000') # '#ffffff' = white

# "blue"   "yellow" "black"  "green"  "brown"  "pink"   "red" 
# interesting_modules_cols = c('#4363d8','#F6D860','#212121','#4E9F3D','#B05E27','pink','#e6194B')

## filter and sort by yael's list
yael = read.table('rnaSeq_main/modules_enrichment/data/SOURCE_Israel_TI_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default/module_function_top10_manually_selected_yael.txt', header = T, sep = '\t')
yael$Yael_short_Numberr[yael$Category == 'ToppCell Atlas'] = yael$Yael_Short

df = df[df$ID %in% yael$ID,]
df = merge(df, yael[,c('Yael_Short','Yael_short_Numberr','ID')], all.x = T, by = 'ID')
df$Name = factor( df$Name, levels = yael$Name )
df$Yael_short_Numberr = factor( df$Yael_short_Numberr, levels = rev(yael$Yael_short_Numberr) )

df$module = factor( df$module, levels = wanted_modules )

# "magenta"
g = ggplot(df) + 
  geom_point(aes(x=module, y = Yael_short_Numberr, colour = Category), shape=15, size=0.1, alpha=0) +
  geom_point(shape=21, aes(x=module, y = Yael_short_Numberr, fill = module, size=FDR_BH)) +
  theme_bw() + ylab('') + 
  # facet_grid(Database~., scales='free') + 
  scale_fill_manual(breaks = wanted_modules, values = wanted_modules)+ 
  theme(axis.ticks.y = element_blank()) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  guides(fill="none", color = guide_legend(override.aes = list(size = 5, alpha=1) ) ) +
  labs(size='-log10(FDR)') + xlab('Module') + 
  scale_colour_brewer(palette = 'Set3') #
# scale_colour_manual(values = unique_cols) #

# Build a legend "bar"
leg <- ggplot(df, aes(y = Yael_short_Numberr, x = 0)) + 
  geom_point(aes(color = Category), shape = 15, size = 3, show.legend = F) + 
  theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm")) +
  # scale_colour_manual(values = unique_cols) #
  scale_colour_brewer(palette = 'Set3') #
# annotate hm with the bar
g = g + annotation_custom(ggplotGrob(leg), 
                          xmin = -.1, xmax = .3, 
                          ymin = 0, ymax =length(unique(df$Name))+0.5) # n + .5 since tiles are 1 in width centered on the value
g


out_path = 'rnaSeq_main/modules_enrichment/res/'
dir.create(out_path)
# ggsave(sprintf('%s/module_function_top%s.pdf', out_path, top_num),plot = g, device = 'pdf', width = 15,height = 30)
ggsave(sprintf('%s/module_function_top%s_yaelFlt.pdf', out_path, top_num),plot = g, device = 'pdf', width = 8,height = 8)

# write.table(x = temp, file = sprintf('%s/module_function_top%s.txt', out_path,top_num), quote = F, sep = '\t',col.names = T, row.names = F)











# df$FDR_BH = -1*log10(df$FDR_BH)
# 
# temp  = df[,c('Name','FDR_BH','module')]
# temp = reshape2::dcast(temp, Name ~ module, value.var = 'FDR_BH')
# row.names(temp) = temp$Name
# temp = temp[,-1]
# temp[is.na(temp)]=0
# 
# # ord = hclust( dist(temp, method = "euclidean") )$order
# # df$Name = factor( df$Name, levels = unique(df$Name)[ord] )
# 
# # ord = hclust( dist(t(temp), method = "euclidean") )$order
# # df$module = factor( df$module, levels = unique(df$module)[ord] )
# g
# df$module = factor(df$module, levels = c('blue','yellow','black','green', 'brown','pink', 'red'))
# # "magenta"
# ggplot(df, aes(x=module, y = Name, fill = module, size=FDR_BH)) + geom_point(shape=21) +
#   theme_bw() + 
#   # facet_grid(Database~., scales='free') + 
#   scale_fill_manual(breaks = interesting_modules, values = interesting_modules)
