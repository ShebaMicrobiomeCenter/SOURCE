source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
library(ggplot2)

source("/pita/users/tzipi/projects/multiomics/SOURCE/multi_analysis/WGCNA/rnaSeq_main/israel_china_met_cmpr.R")
col2 = colorRampPalette(rev( c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                               "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061") ))(200)

mdc = mets_all
mdc$variable = make.names(mdc$variable)
mdc = mdc[mdc$variable %in% mic$variable[mic$both_cmp == 'Same ditrection'],]

wanted_modules = c('yellow','green','red','pink','purple','tan','salmon','black','brown')

mdc$module = factor(mdc$module, levels = wanted_modules)

temp = mdc[mdc$country == 'Israel' & mdc$module == 'yellow',]
mdc$variable = factor(mdc$variable, levels = temp$variable[order(temp$correlation)]   )
mdc$country = factor(mdc$country, levels = c('Israel','China'))

heatmap_ic = ggplot(mdc, aes(x=module,y=variable,fill = correlation)) +
  geom_tile() + xlab('Module') +
  scale_fill_gradientn(colours = col2, limits = c(-1,1)) +
  # geom_point(aes(shape = p_value<=0.05)) + scale_shape_manual(values = c(NA, 8)) +
  # geom_text(aes(label=sprintf('%.2f\n%.2e',correlation,p_value )), size=3) +
  # geom_text(aes(label=sprintf('%.2e',p_value )), size=3) +
  # geom_point(aes(alpha = q_value <= 0.25), size=1, shape = 16) +
  # scale_alpha_manual(values = c(0,1)) + 
  # scale_colour_manual(values = c('black','white')) +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size= 10)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size= 10)) + coord_flip() +
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  facet_grid(~country, scales = 'free') + 
  # ggtitle('China metabolomics all') +
  # theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  ylab('Metabolite')
ggsave('rnaSeq_main/final_figs/israle_china_shared_dir_mets_heatmap.pdf',plot = heatmap_ic, device = 'pdf', width = 7,height = 4.5)

# heatmap_ic = heatmap_ic + geom_point(aes(alpha = q_value <= 0.25), size=1, shape = 16) +
#   scale_alpha_manual(values = c(0,1)) 
# ggsave('rnaSeq_main/final_figs/israle_china_shared_dir_mets_heatmap_q.pdf',plot = heatmap_ic, device = 'pdf', width = 7,height = 4.5)

heatmap_ic_small = heatmap_ic + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  
ggsave('rnaSeq_main/final_figs/israle_china_shared_dir_mets_heatmap_small.pdf',plot = heatmap_ic_small, device = 'pdf', width = 7,height = 2.2)
