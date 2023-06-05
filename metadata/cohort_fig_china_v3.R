library(ggplot2)

## China

## rnaseq
rnaSeq_metadata_file = '../../rnaSeq/China_Israel/data/SOURCE_israel_china_rnaSeq_map_v2.txt'
rna_map = read.table(file = rnaSeq_metadata_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'))
rna_map = rna_map[rna_map$Country == 'China',]

name = 'SOURCE_China_TI_ageMatch_withFFQv10_oneSample_pwr_14_netType_signed_hybrid_minMdSz_30_WGCNA_default'
ti = load( sprintf('../../multi_analysis/WGCNA/rnaSeq_main/res_v2/%s/data.RData', name) )
rna_map_ti = rna_map[make.names(rna_map$SampleID) %in% make.names(row.names(metadata_df)),]
name = 'SOURCE_China_R_ageMatch_withFFQv10_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default'
israel_r = load( sprintf('../../multi_analysis/WGCNA/rnaSeq_main/res_v2/%s/data.RData', name) )
rna_map_r = rna_map[make.names(rna_map$SampleID) %in% make.names(row.names(metadata_df)),]
rna_map = rbind(rna_map_ti, rna_map_r)

## ffq
ffq_file = '../../metadata/FFQ/FFQ_v10.txt'
ffq = read.table(ffq_file, header = T, sep = '\t')
names(ffq)[1] = 'pn_ID'
ffq = ffq[!grepl('A',ffq$pn_ID),]

## 16s
taxa_map_file = '/pita/users/tzipi/projects/multiomics/SOURCE/16S/data/SOURCE_16S_china_map_v4_Stool.txt'
taxa_map = read.table(taxa_map_file, header = T, sep = '\t')

## mtg
mtg_map_file = '../../metagenomics/data/source_chinese_mtg_species.tsv'
mtg_map = read.table(mtg_map_file, header = T, sep = '\t')

## metabolomics
mtb_map_file = '../../metabolomics/China/data/SOURCE_china_mets_table_norm_v2.tsv'
mtb_map = read.table(mtb_map_file, header = T, sep = '\t', quote = '', comment.char = '')


## metadata
metadata_map_file = '/pita/users/tzipi/projects/multiomics/SOURCE/metadata/alona_data/SOURCE_Israel_China_data_v9.txt'
dx_map = read.table(metadata_map_file, header = T, sep = '\t')

pns = unique( c(rna_map$pn_ID, ffq$pn_ID, taxa_map$pn_ID) )

df = data.frame(pn_ID = pns)
# df = merge(df, dx_map[,c('pn_ID','Patient_group')], by='pn_ID',all.x = T)
df = merge(df, dx_map[,c('pn_ID','Patient_group2')], by='pn_ID',all.x = T); names(df)[2] = 'Patient_group'

df$rnaSeq_TI = df$pn_ID %in% rna_map$pn_ID[rna_map$Location == 'TI'] 

df$ffq = df$pn_ID %in% ffq$pn_ID

df$taxa_Stool = df$pn_ID %in% taxa_map$pn_ID[taxa_map$Location_inflamation == 'Stool'] 
df$Metagenomics = df$pn_ID %in% names(mtg_map)

df$Metabolomics = df$pn_ID %in% names(mtb_map)

df$Demographics = df$pn_ID %in% dx_map$pn_ID[!is.na(dx_map$Age_years) | !is.na(dx_map$Gender)] 
df$Environment = df$pn_ID %in% dx_map$pn_ID[dx_map$IOIBD_completed == 'Yes'] 
df$CRP = df$pn_ID %in% dx_map$pn_ID[ !is.na(dx_map$CRP_numeric) ]
# df$Calprotectin = df$pn_ID %in% dx_map$pn_ID[ !is.na(dx_map$Calprotectin_numeric) ] 

df$fill_sum = rowSums(df[,3:dim(df)[2]])

df$Patient_group[df$Patient_group == 'Urban_health'] = 'China urban'
df$Patient_group[df$Patient_group == 'Rural_health'] = 'China rural'
df$Patient_group[df$Patient_group == 'Chinese_crohns'] = 'China CD'
df$Patient_group[df$Patient_group == 'Rural_health_<50%_in_city'] = 'China rural'
df$Patient_group[df$Patient_group == 'Rural_health_>50%_in_city'] = 'China rural-urban'


dfm = reshape::melt(data = df, id.vars = c('pn_ID','fill_sum','Patient_group'), variable_name = 'Data_type')

dfm$Data_type = gsub('taxa','16S',dfm$Data_type)
dfm$Data_type = gsub('ffq','FFQ',dfm$Data_type)

dfm$pn_ID = factor(dfm$pn_ID, levels = df$pn_ID[order(df$fill_sum)])
dfm$Type = gsub('_.*','',dfm$Data_type)
dfm$Source = gsub('.*_','',dfm$Data_type)

dfm$Type[dfm$Type %in% c('Demographics','Environment','CRP','Calprotectin')] = 'Metadata'
dfm$Type[dfm$Type %in% c('16S','Metagenomics')] = 'Bacteria'
dfm$Type[dfm$Type %in% c('rnaSeq')] = 'RNA-Seq'
dfm$Type[dfm$Type %in% c('metabolomics')] = 'Metabolomics'

dfm$Type[dfm$Type %in% c('Demographics','Environment','CRP','Calprotectin')] = 'Metadata'
dfm$Type[dfm$Type %in% c('16S','Metagenomics')] = 'Bacteria'

dfm$Data_type[dfm$Data_type %in% c('rnaSeq_TI')] = 'RNA-Seq TI'
dfm$Data_type[dfm$Data_type %in% c('Metagenomics')] = 'Metagenomics stool'
dfm$Data_type[dfm$Data_type %in% c('Metabolomics')] = 'Metabolomics stool'
dfm$Data_type = gsub('_',' ',dfm$Data_type)

dfm$Type = factor(dfm$Type, levels = c('Metadata','FFQ','Bacteria','RNA-Seq','Metabolomics'))
dfm$Data_type = factor(dfm$Data_type, levels = c('Demographics','Environment','CRP','Calprotectin','FFQ','16S Stool','Metagenomics stool','RNA-Seq TI','Metabolomics stool'))

dfm$Patient_group3 = dfm$Patient_group
for (pg in unique(dfm$Patient_group))
{
  dfm$Patient_group3[dfm$Patient_group == pg] = sprintf('%s\nn=%s', gsub('_',' ',pg), length(unique( dfm$pn_ID[dfm$Patient_group == pg] )))
}

# dfm$Patient_group3 = factor(dfm$Patient_group3, levels = c('China CD\nn=40','China urban\nn=121','China rural\nn=162'))
dfm$Patient_group3 = factor(dfm$Patient_group3, levels = c('China CD\nn=40','China urban\nn=121','China rural-urban\nn=74','China rural\nn=88'))
g = ggplot(dfm, aes(x=Data_type, y=pn_ID, alpha = value, fill = Type)) + 
  geom_tile() + 
  scale_alpha_manual(values = c(.9,0), name = 'Data\navailable', breaks = c(T,F)) +
  # scale_fill_manual(values = c('#006E7F','#F8CB2E','#EE5007','#990000')) + 
  scale_fill_manual(values = c('#293462','#006E7F','#F8CB2E','#EE5007','#990000')) + 
  # ggtitle(sprintf('additional %s FFQ and %s 16S stool samples\nChina',sum(df_missing$ffq), sum(df_missing$taxa_Stool))) + 
  ggtitle('China') + 
  ylab('Subject') + xlab('') + 
  # facet_grid(~Type + Source, scales = 'free') + 
  facet_grid(Patient_group3~. , scales = 'free', space = 'free') + 
  theme_bw() +scale_y_discrete(expand=c(0,0)) +scale_x_discrete(expand=c(0,0)) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #, axis.text.y = element_blank()
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.key = element_rect(colour = "black")) 


ggsave('China_cohort_figure_full_g3.pdf',plot = g, device = 'pdf', width = 4,height = 35, limitsize = FALSE)

g2 = g+theme(axis.text.y = element_blank(), axis.ticks.y = element_blank() )
ggsave('China_cohort_figure_g3.pdf',plot = g2, device = 'pdf', width = 3.5,height = 10, limitsize = FALSE)

row.names(dx_map) = dx_map$pn_ID
dx_map = dx_map[df$pn_ID,]
df$Age = dx_map$Age_years
df$Gender = dx_map$Gender
df$BMI = dx_map$BMI
df$CRP_mg_L = dx_map$CRP_numeric
df$Calprotectin_ug_g = dx_map$Calprotectin_numeric

write.table(x = df, file = 'China_basic_map.txt', quote = F, sep = '\t',row.names = F)


