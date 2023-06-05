library(WGCNA)
library(ggplot2)
options(stringsAsFactors = FALSE);

source('/pita/users/tzipi/code/R_figs/WGCNA_funcs.R')

min_p = 0.25
data_type = 'all'

# cor_type = 'spearman'
cor_type = 'WGCNA_default'

in_path = 'rnaSeq_main/res_v2/'
name = sprintf('SOURCE_Israel_TI_withStoolMetablolomicsV2LogClean_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_%s',cor_type)
israel_ti = load( sprintf('%s/%s/data.RData',in_path, name) )



metadata_prms = c('Gender','Age_years','BMI',
                  'Active_Smoker_manual','Previous_or_Current_Tobacco_Use_manual',
                  'Disease','inflammed','CRP_numeric','Calprotectin_numeric',# 'CRP_log','Calprotectin_log',
                  'health_index_same_sample','health_index_stool')

## heatmap using predefined module list
wanted_modules = c('yellow','green','red','pink', 'purple','tan','salmon','black','brown')
wanted_modules = sprintf('ME%s',wanted_modules)

res = WGCNA_ME_met_heatmap_full_res(moduleTraitCor[wanted_modules,], moduleTraitPvalue[wanted_modules,], name, metadata_prms, min_p = min_p, wanted_modules_order = wanted_modules, fdr_flag_non_metadata = T, show_p_value_fdr_flag = T )
mdc = res[[3]]
mdc$variable = make.names(mdc$variable)

col2 = colorRampPalette(rev( c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                               "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061") ))(200)
mdc = res[[2]]
mdc = mdc[!mdc$variable %in% metadata_prms,]

mdc$Country = 'Israel'
mdc$var_type = 'Metabolites'

## changing clustering, unnecessary here.
# dfm = as.data.frame(maditr::dcast(data = mdc,formula = variable~module, 
#                                   fun.aggregate = sum,value.var = "correlation") )
# row.names(dfm) = dfm$variable; dfm = dfm[,-1]
# ht = heatmap(as.matrix(dfm, scale='none', col=col2, symm = F, scale='none'))
# # mdc$variable = factor(mdc$variable, levels = row.names(dfm)[ht$rowInd])

g_full = ggplot(mdc, aes(x=module,y=variable,fill = correlation)) +
  geom_tile() + xlab('Module') + 
  scale_fill_gradientn(colours = col2, limits = c(-1,1)) +
  # geom_point(aes(shape = p_value<=0.05)) + scale_shape_manual(values = c(NA, 8)) +
  # geom_text(aes(label=sprintf('%.2f\n%.2e',correlation,p_value )), size=3) + 
  scale_colour_manual(values = c('black','white')) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size= 10)) + 
  # theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size= 10)) + coord_flip() +
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) + ggtitle('All') + 
  ylab('Metabolite')
g = g_full + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) 
g_full = g_full + geom_text(aes(label=signif(p_value, 1), colour =p_value < 0.001 ), size=3) 
  
# g = g + geom_point(aes(x=module,y=variable, alpha = q_value <= 0.25), shape = 20, size = 0.1) + 
#   scale_alpha_manual(values = c(0,1), name= 'Q <= 0.25')#   + 
#   # facet_grid(var_type + Country~., scales = 'free', space = 'free')

out_path = 'rnaSeq_main/final_figs/'
ggsave(sprintf('%s/module_cor_heatmap_israel_mets_q%s_small_withSig.pdf',out_path, min_p),plot = g, device = 'pdf', width = 6,height = 3, limitsize = F)


y_lables = ggplot_build(g)$layout$panel_params[[1]]$y$get_labels()

## for looking at only a subcohort ( only diease/healthy)
subcohort_name = 'Disease'
wanted_pos = metadata_df$Disease == 1

res = WGCNA_calculate_module_eigengenes_and_corrs(ftr_df[wanted_pos,],moduleColors, clean_metadata_eampty_noVar( metadata_df[wanted_pos,] ), cor_type, name)
moduleTraitCor_d = res[[2]]
moduleTraitPvalue_d = res[[3]]
wanted_labels = y_lables[y_lables %in% colnames(moduleTraitCor_d)]

res = WGCNA_ME_met_heatmap_full_res(moduleTraitCor_d[wanted_modules,wanted_labels], moduleTraitPvalue_d[wanted_modules,wanted_labels], name, metadata_prms, wanted_modules_order = wanted_modules, cluster_cols_order_flag = F )

mdc2 = res[[2]]
mdc2 = mdc2[!mdc2$variable %in% metadata_prms,]
mdc2$q_value = p.adjust(mdc2$p_value, method = 'BH', n = length(mdc2$p_value))
g2_full = ggplot(mdc2, aes(x=module,y=variable,fill = correlation)) +
  geom_tile() + xlab('Module') + 
  scale_fill_gradientn(colours = col2, limits = c(-1,1)) +
  # geom_point(aes(shape = p_value<=0.05)) + scale_shape_manual(values = c(NA, 8)) +
  # geom_text(aes(label=sprintf('%.2f\n%.2e',correlation,p_value )), size=3) + 
  scale_colour_manual(values = c('black','white')) + 
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size= 10)) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size= 10)) + coord_flip() +
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) + 
  ggtitle('Disease') + 
   ylab('Metabolite')
g2 = g2_full + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) 
g2_full = g2_full + geom_text(aes(label=signif(p_value, 1), colour =p_value < 0.001 ), size=3) 
# g2 = g2 + geom_point(aes(x=module,y=variable, alpha = q_value <= 0.25), shape = 20, size = 0.1) + 
#   scale_alpha_manual(values = c(0,1), name= 'Q <= 0.25')#   + 

library(patchwork)
g3 = (g + theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust=1)))+
  (g2 + ylab(''))

# ggsave(sprintf('%s/module_cor_heatmap_israel_mets_q%s_small_withSig_all_disease.pdf',out_path, min_p),plot = g3, device = 'pdf', width = 5,height = 5.5, limitsize = F)
ggsave(sprintf('%s/module_cor_heatmap_israel_mets_q%s_smallall_disease.pdf',out_path, min_p),plot = g3, device = 'pdf', width = 5,height = 5.5, limitsize = F)

p  = ( g_full + theme(legend.position = 'none') + 
         ggtitle('Israel all CD and controls') ) + 
  ( g2_full + theme( axis.text.y = element_blank(), axis.ticks.y = element_blank() ) + 
      ggtitle('Israel only within CD') + guides(color = "none") + ylab('') )
ggsave(sprintf('%s/module_cor_heatmap_israel_all_and_CD_allq%s.pdf',out_path, min_q),plot = p, device = 'pdf', width = 12,height = 25)


mdc$group = 'CD and control'
mdc2$group = 'CD'
mdc2$Country = 'Israel'
mdc2$var_type = 'Metabolites'
temp = rbind(mdc, mdc2)
write.table(x = temp, 
            file = sprintf('%s/module_cor_heatmap_israel_mets_q%s.txt',out_path , min_p),
            quote = F, sep = '\t', row.names = F)

mdc$id = sprintf('%s_%s',mdc$module, mdc$variable)
mdc2$id = sprintf('%s_%s',mdc2$module, mdc2$variable)
mdc = mdc[order(mdc$id),]
mdc2 = mdc2[order(mdc2$id),]

t = mdc
t$cd_cor = mdc2$correlation
t=t[t$q_value<=0.25,]
t$same_dir = ifelse( (t$correlation < 0 & t$cd_cor < 0) | (t$correlation > 0 & t$cd_cor > 0),T,F )
print( sum(t$same_dir)/length(t$same_dir) )
print(table(t$same_dir))

print(table(t$cd_cor>0))
p_cd = sum( t$cd_cor>0) / length(t$cd_cor) 
print(table(t$correlation>0))
p_all = sum( t$correlation>0) / length(t$correlation) 
prob = p_cd*p_cd + (1-p_cd)*(1-p_all)
binom.test(x = sum(t$same_dir==T), n = length(t$cd_cor) , p = prob)



