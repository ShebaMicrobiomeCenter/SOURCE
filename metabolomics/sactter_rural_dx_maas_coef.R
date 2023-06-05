source('/pita/users/tzipi/code/R_figs/maaslin2_funs.R')
source('/pita/users/tzipi/code/R_figs/figs_funcs.R')
library(Maaslin2)
library(ggplot2)

taxonomy_flag = F
# taxonomy_file = '/pita/users/tzipi/projects/16S/all_merges/DB1-20_merge/res/16S_DB1-20_merged/taxonomy/metadata_v2_source.tsv'

name = 'source'
path = '.'
taxa_file = '../data/SOURCE_china_mets_table_norm_v2.tsv'
na_str = c('no_data','_','NA','unknown', 'other','na','No_followup')
map = read.table(file = '../data/SOURCE_china_mets_map.tsv', header = T,sep='\t', quote = '',comment.char = '')
taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_str, row.names = 1)
taxa = as.data.frame(t(taxa))
names(taxa) = make.names(names(taxa))

sig_var = 'Patient_group2'
check_val = sig_var

head_path = sprintf('%s/res/maaslin2_v2/',path)

name = 'gender_age_group_inRural'
maas_path = sprintf('%s/%s',head_path, name)
mas_res_r = read.table(sprintf('%s/res/all_results.tsv',maas_path), header = T, quote = '',comment.char = '')

name = 'gender_age_group_inUrban'
maas_path = sprintf('%s/%s',head_path, name)
mas_res_d = read.table(sprintf('%s/res/all_results.tsv',maas_path), header = T, quote = '',comment.char = '')

qval_cut = 0.25

mas_res_r = mas_res_r[mas_res_r$metadata == sig_var & mas_res_r$qval<=0.25,]
mas_res_d = mas_res_d[mas_res_d$metadata == sig_var & mas_res_d$feature %in% mas_res_r$feature,]
mas_res_r = mas_res_r[order(mas_res_r$feature),]
mas_res_d = mas_res_d[order(mas_res_d$feature),]

df = data.frame(Feature = mas_res_r$feature, 
                coef_rural = mas_res_r$coef, coef_dx = -1*mas_res_d$coef,
                q_rural = mas_res_r$qval, q_dx = mas_res_d$qval)

# df$FC_rural = 0
# df$FC_dx = 0
# for ( i in 1:length(df$Feature) )
# {
#   ru_mean = taxa[ map$SampleID[map$Patient_group2 == 'Rural_health_>50%_in_city'], df$Feature[i] ]
#   r_mean = taxa[ map$SampleID[map$Patient_group2 == 'Rural_health_<50%_in_city'], df$Feature[i] ]
#   df$FC_rural[i] = mean(ru_mean) / mean(r_mean)
#   u_mean = taxa[ map$SampleID[map$Patient_group2 == 'Urban_health'], df$Feature[i] ]
#   cd_mean = taxa[ map$SampleID[map$Patient_group2 == 'Chinese_crohns'], df$Feature[i] ]
#   df$FC_dx[i] = mean(u_mean) / mean(cd_mean)
# }

Direction = ifelse(df$coef_dx > 0,'Higher in\nCD\nRural-urban\n','Higher in\nUrban controls\nRural\n')
Direction[df$q_dx>0.25] = 'Not significant\nin CD vs\nRural controls\n'
Feature2 = df$Feature
# remove long names
Feature2[Feature2 %in% c(# 'Ketoleucine',
                         'cis.7.10.13.16.Docosatetraenoic.acid',
                         'Coumaric.acid.4.Hydroxycinnamic.acid',
                         'X3.Methyl.2.oxopentanoic.acid',
                         'X2.3.Diaminopropionic.acid',
                         'Pyrrole.2.carboxylic.acid',
                         'X10Z.Nonadecenoic.acid')] = NA
library(ggrepel)
g = ggplot( df, 
            # aes(x=coef_rural, y=coef_dx, colour = q_dx<=0.25) ) + 
            aes(x=coef_rural, y=coef_dx, label = Feature2) ) +
  geom_vline(xintercept = 0, colour ='gray80', linetype="dashed") + 
  geom_hline(yintercept = 0, colour ='gray80', linetype="dashed") + 
  scale_colour_manual(values = c('red','blue','gray60'), name = '') +
  # geom_text_repel(size = 2, seed = 4, nudge_y = 0.2, box.padding = 7) + 
  geom_point(aes(colour = Direction)) + theme_bw() + theme(panel.grid = element_blank()) + 
  # geom_text_repel(size = 2, seed=2023,nudge_y = 0.3 ,point.padding = 11) +
            # aes(x=FC_rural, y=FC_dx, colour = q_dx<=0.25) ) + 
  xlab('Rural-urabn vs rural coef') + ylab('CD vs urban controls coef') 
ggsave('res/rural_dx_maas_overlap_scatter.tiff',
       plot = g, device = 'tiff', width =4,height = 2.5, compression='lzw', points= 0.8)

cor.test(df$coef_rural, df$coef_dx, test_type='spearman')

# taxa = taxa[taxa$Dx == 'healthy',]

# mas_taxa = read.table(file = sprintf('%s/taxonomy.tsv',maas_path), header = T,sep='\t', row.names = 1, quote = '',comment.char = '')
# # mas_met = read.table(file = sprintf('%s/metadata.tsv',maas_path), header = T,sep='\t', quote = '',comment.char = '')
# mas_met = read.table(file = '../data/SOURCE_china_mets_map.tsv', header = T,sep='\t', quote = '',comment.char = '')

# row.names(taxa) = taxa$SampleID
# taxa = taxa[row.names(mas_taxa),]
# taxa = mas_taxa

temp = mas_res_u


mas_res_r_f = mas_res_r
mas_res_r_f = mas_res_r_f[mas_res_r_f$metadata == sig_var,]
mas_res_u_f = mas_res_u
mas_res_u_f = mas_res_u_f[mas_res_u_f$metadata == sig_var,]

mas_res_r_f$type ='Rural'
mas_res_u_f$type ='Dx'

mas_res = rbind(mas_res_r_f, mas_res_u_f)
# mas_res = mas_res[mas_res$qval <= 0.25, ]
mas_res = mas_res[mas_res$feature %in% 
                    mas_res$feature[ mas_res$type == 'Rural' & mas_res$qval <= 0.25], ]



# out_path = sprintf('%s/res/my_figures/', maas_path)
# dir.create(out_path)
# # out_path = 'art_plots/res/'


# library(ggrepel)
# 
# col = ifelse(mas_res_f$qval <= qval_cut, 'Up','NS')
# col[col == 'Up' & mas_res_f$coef<0] = 'Down'
# 
# lbl = gsub('.*__','',mas_res_f$feature)
# lbl[col == 'NS'] = NA
# lbl = gsub('_','\n',lbl)
# volc_p = ggplot(mas_res_f, aes(x=coef, y=-log10(qval), label = lbl)) + 
#   geom_point(aes(colour = col)) + 
#   # geom_text_repel(size = 2, seed = 4, nudge_y = 0.2, max.overlaps = 2, box.padding = 0.1) + 
#   geom_text_repel(size = 2, seed = 4, nudge_y = 0.2, box.padding = 5) + 
#   scale_color_manual(values = c('blue','gray','red'), name='') + 
#   theme_bw() + theme(panel.grid = element_blank())
# ggsave(sprintf('%s/%s_q%s_volcano.tiff', out_path, check_val, qval_cut),
#        plot = volc_p, device = 'tiff', width =4,height = 2.5, compression='lzw', points= 0.8)
# 
# 
# # df = mas_taxa
# df = taxa
# 
# df = df[,mas_res_f$feature[ mas_res_f$qval <= qval_cut ]]
# df$SampleID = row.names(df)
# dfm = melt(df, id.vars = 'SampleID', variable_name = 'feature')
# dfm = merge(x = dfm, y = mas_met[,c('SampleID',sig_var)], by = 'SampleID')
# # dfm = merge(x = dfm, y = mas_res_f[,c('feature','coef','pval','qval')], by = 'feature')
# dfm$value = dfm$value + 0.0001
# 
# mas_res_f2 = mas_res_f
# mas_res_f2 = mas_res_f2[mas_res_f2$feature %in% dfm$feature, ]
# mas_res_f2 = mas_res_f2[order(mas_res_f2$coef),]
# dfm$feature = factor(dfm$feature, levels = mas_res_f2$feature)
# 
# mas_res_f2$ord = vector(mode = 'numeric',length = length(mas_res_f2$feature))
# for ( i in 1:length(mas_res_f2$feature) )
# {
#   mas_res_f2$ord[i] = mean( dfm$value[dfm$feature == mas_res_f2$feature[i]  & 
#                                         dfm$Patient_group2 %in% c('Rural_health_<50%_in_city',
#                                                                   'Rural_health_>50%_in_city')]  )
#   if ( mas_res_f2$coef[i] < 0 )
#     mas_res_f2$ord[i] = -1 * mas_res_f2$ord[i]
# }
# dfm$feature = factor(dfm$feature, levels = mas_res_f2$feature[order( mas_res_f2$ord )])
# dfm$Patient_group2 = factor(dfm$Patient_group2, levels =
#                               c('Rural_health_<50%_in_city','Rural_health_>50%_in_city',
#                                 'Urban_health','Chinese_crohns'))
# dfm$group = ifelse(dfm$Patient_group2 %in% c('Urban_health','Chinese_crohns'),'Urban','Rural')
# 
# box_p = ggplot(dfm) + 
#   geom_boxplot(aes(x=feature, y= value, fill =  .data[[sig_var]])) +
#   # scale_fill_manual(values = c('#7FB77E','#F7A76C'), name = '',
#   #                   labels = c('Rural_health_<50%_in_city' = 'Rural',
#   #                              'Rural_health_>50%_in_city' = 'Rural-Urban')) +
#   scale_fill_manual(values = c('#7FB77E','#F7A76C','#3DB2FF','#EB5353'), name = '',
#                     labels = c('Rural_health_<50%_in_city' = 'Rural',
#                                'Rural_health_>50%_in_city' = 'Rural-Urban',
#                                'Urban_health' = 'Urban',
#                                'Chinese_crohns' = 'CD')) +
#   # scale_fill_manual(values = c('#82CD47','#1F8A70'), name = '') +
#   xlab('') + ylab('Relative Abundance (RA)') + 
#   scale_y_log10() + theme_bw() + coord_flip() +
#   facet_grid(~group, scales= 'free') + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7),
#         panel.grid.minor = element_blank())
# # theme(axis.text.x = element_text(angle = 45, hjust=1))
# 
# box_p = format_fig(fig = box_p, set_ttl_flag = F) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7),
#         panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank()) 
# # ggsave(sprintf('%s/%s_q%s_boxplot.tiff', out_path, check_val, qval_cut),plot = box_p, device = 'tiff', width =9,height = 5, compression='lzw')
# # ggsave(sprintf('%s/%s_q%s_boxplot.tiff', out_path, check_val, qval_cut),
# #        # plot = box_p, device = 'tiff', width =8,height = 6, compression='lzw')
# #        plot = box_p, device = 'tiff', width =8,height = 12, compression='lzw')
