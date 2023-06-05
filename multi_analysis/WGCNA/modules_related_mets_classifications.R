source('/pita/users/tzipi/code/R_figs/figs_funcs.R')

mets_types_file = '../../metabolomics/run2_oct22/data/metabolites_hmdb_classification_nina.txt'
mets_types = read.table(mets_types_file, header = T, sep = '\t', quote = '', comment.char = '')
row.names(mets_types) = make.names(mets_types$Feature)

wanted_modules = c('yellow','green','red','pink', 'purple','tan','salmon','black','brown')
mets_path = 'rnaSeq_main/module_related_mets/'

bad_mets = c()
all_df = data.frame(met = c(), module = c(), direction = c(), final_class = c())
for (md in wanted_modules)
{
  for (dir in c('up','down'))
  {
    met_file = sprintf('%s/%s_%s.txt', mets_path, md,dir)
    if ( file.size(met_file) != 0L ) 
    {
      m = read.table(met_file, header = F, sep = '\t', quote = '', comment.char = '')
      print(sprintf('%s %s', md, dir))
      print(table(m$V1 %in% make.names(mets_types$Feature)))
      if ( ! all(m$V1 %in% make.names(mets_types$Feature)) )
      {
        print(m$V1[ ! m$V1 %in% make.names(mets_types$Feature) ])
        bad_mets = c(bad_mets,m$V1[ ! m$V1 %in% make.names(mets_types$Feature) ] )
      }
      good_mets = m$V1[ m$V1 %in% make.names(mets_types$Feature) ]
      df = data.frame(met = good_mets, module = md, direction = dir, 
                      final_class =  mets_types[good_mets,]$Final_Class)
      all_df = rbind(all_df, df)
      
    }
  }
}
bad_mets = unique(bad_mets)

all_df$module = factor(all_df$module, levels = wanted_modules)
all_df$sum = 1

disease_modules = c('purple','tan','salmon','black','brown')
all_df$disease_direction = all_df$direction
all_df$disease_direction[all_df$module %in% disease_modules & all_df$direction == 'up'] = 'down'
all_df$disease_direction[all_df$module %in% disease_modules & all_df$direction == 'down'] = 'up'
all_df$disease_direction = ifelse(all_df$disease_direction=='up','Control associated','Disease associated')


cols = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
gm= ggplot(data = all_df, aes(x=module, y=sum, fill = final_class)) + 
  geom_bar(stat="identity", alpha = 0.8) + 
  facet_grid(disease_direction~.) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid = element_blank()) + 
  scale_fill_manual(values = cols)
ggsave(filename = 'rnaSeq_main/module_related_mets/modules_related_mets_classification_bar_by_module.pdf',
       plot = gm, device = 'pdf',dpi = 600, height = 5, width = 7)

wanted_modules = c('yellow','green','red',
                   'pink','purple','tan',
                   'salmon','black','brown')
wanted_modules_cols = c('#F7C04A','#539165','#EB455F',
                        '#F7C8E0','#645CBB','#d2b48c',
                        '#FA8072','#2C3333','#A86464')
g = ggplot(data = all_df, aes(x=final_class, y=sum, fill = module)) + 
  geom_bar(stat="identity") + 
  # facet_grid(direction~., scales = 'free') + 
  facet_grid(disease_direction~.) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid = element_blank()) + 
  scale_fill_manual(values = wanted_modules_cols)
ggsave(filename = 'rnaSeq_main/module_related_mets/modules_related_mets_classification_bar.pdf',
       plot = g, device = 'pdf',dpi = 600, height = 6, width = 5)
