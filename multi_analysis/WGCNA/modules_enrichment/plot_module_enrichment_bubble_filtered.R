library(ggplot2)

wanted_modules = c('yellow','green','red',
                   'pink','purple','tan',
                   'salmon','black','brown')
wanted_modules_cols = c('#F7C04A','#539165','#EB455F',
                        '#F7C8E0','#645CBB','#d2b48c',
                        '#FA8072','#2C3333','#A86464')
top_num = 10

df = read.table('rnaSeq_main/modules_enrichment/data/SOURCE_Israel_TI_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default/module_function_top10_manually_selected_yael_mar23.txt', header = T, sep = '\t')
df$Yael_short_Numberr[df$Category == 'ToppCell Atlas'] = df$Yael_Short[df$Category == 'ToppCell Atlas']

df$Yael_short_Numberr = gsub('_',' ',df$Yael_short_Numberr)

df = df[!is.na(df$include_YES),c('new_order','Yael_short_Numberr','Category','ID',
                              wanted_modules)]

df_old = df
df = reshape::melt(data = df, 
                   id.vars =  c('new_order','Yael_short_Numberr','Category','ID'),
                   variable_name = 'module')
names(df)[names(df) == 'value']  = 'FDR_BH'
df$FDR_BH[df$FDR_BH <=0.05] = NA

df$module = factor( df$module, levels = wanted_modules )
df$Yael_short_Numberr = factor( df$Yael_short_Numberr, levels = rev(df_old$Yael_short_Numberr) )


df$Category = factor(df$Category,
                         levels = c('GO: Biological Process','GO: Cellular Component',
                                    'GO: Molecular Function','Pathway',
                                    'ToppCell Atlas', 'Disease'))
g = ggplot(df) + 
  geom_point(aes(x=module, y = Yael_short_Numberr, colour = Category), shape=15, size=0.1, alpha=0) +
  geom_point(shape=21, aes(x=module, y = Yael_short_Numberr, fill = module, size=FDR_BH)) +
  theme_bw() + ylab('') + 
  # facet_grid(Database~., scales='free') + 
  scale_fill_manual(breaks = wanted_modules, values = wanted_modules_cols)+ 
  theme(axis.ticks.y = element_blank(), axis.text.y = element_text(size=8),
        panel.grid.major.x = element_blank()) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  guides(fill="none", color = guide_legend(override.aes = list(size = 5, alpha=1) ) ) +
  labs(size='-log10(FDR)') + xlab('Module') + 
  scale_colour_brewer(palette = 'Set3') #
# scale_colour_manual(values = unique_cols) #

df_old$Category = factor(df_old$Category,
                         levels = c('GO: Biological Process','GO: Cellular Component',
                                    'GO: Molecular Function','Pathway',
                                    'ToppCell Atlas', 'Disease'))
df_old$Yael_short_Numberr = factor( df_old$Yael_short_Numberr, levels = rev(df_old$Yael_short_Numberr) )

# Build a legend "bar"
leg <- ggplot(df_old, aes(y = Yael_short_Numberr, x = 0)) + 
  geom_point(aes(color = Category), shape = 15, size = 3, show.legend = F) + 
  theme_classic() + 
  scale_colour_brewer(palette = 'Set3') +
  theme(axis.title = element_blank(), axis.line = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm")) 
  # scale_colour_manual(values = unique_cols) #

# annotate hm with the bar
g = g + annotation_custom(ggplotGrob(leg), 
                          xmin = -.1, xmax = .4, 
                          ymin = 0, ymax =length(unique(df_old$Yael_short_Numberr))+0.5) # n + .5 since tiles are 1 in width centered on the value
g

out_path = 'rnaSeq_main/modules_enrichment/res/'
ggsave(sprintf('%s/module_function_top%s_yaelFlt_mar23.pdf', out_path, top_num),
       plot = g, device = 'pdf', width = 7.7,height = 4.3)


