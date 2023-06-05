library(ggplot2)

## israel

## rnaSeq
rnaSeq_metadata_file = '../../rnaSeq/China_Israel/data/SOURCE_israel_china_rnaSeq_map_v2.txt'
rna_map = read.table(file = rnaSeq_metadata_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'))
rna_map = rna_map[rna_map$Country == 'Israel',]

name = 'SOURCE_Israel_TI_withStoolMetablolomicsV2LogClean_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default'
israel_ti = load( sprintf('../../multi_analysis/WGCNA/rnaSeq_main/res_v2/%s/data.RData', name) )
rna_map_ti = rna_map[make.names(rna_map$SampleID) %in% make.names(row.names(metadata_df)),]
name = 'SOURCE_Israel_R_withFFQv10_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default'
israel_r = load( sprintf('../../multi_analysis/WGCNA/rnaSeq_main/res_v2/%s/data.RData', name) )
rna_map_r = rna_map[make.names(rna_map$SampleID) %in% make.names(row.names(metadata_df)),]
rna_map = rbind(rna_map_ti, rna_map_r)

name = 'SOURCE_Israel_TI_withStoolMetablolomicsV2LogClean_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default'
israel_ti = load( sprintf('../../multi_analysis/WGCNA/rnaSeq_main/res_v2/%s/data.RData', name) )
rna_map_ti = rna_map[make.names(rna_map$SampleID) %in% make.names(row.names(metadata_df)),]

## metabolomics
# met_metadata_file = '/pita/users/tzipi/projects/multiomics/nina_metabolomics/data/nina_10may22/Stool_metadata_v2.txt'
met_metadata_file = '../../metabolomics/run2_oct22/data/SOURCE2022_Stool_metadata.txt'
met_map = read.table(file = met_metadata_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'), comment.char = '', quote="\"")
# met_map = met_map[met_map$Cohort == 'SOURCE',]
# met_map$pn_ID = gsub('S','A',met_map$pn_ID)
## filter to what passed QC
met_file = '../../metabolomics/run2_oct22/data/SOURCE2022_Stool_Data_MZM13_N2_filtered_v2.tsv'
met = read.table(file = met_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'), comment.char = '', quote="\"")
met_map = met_map[make.names(met_map$SampleID) %in% names(met),]

## ffq
ffq_file = '../../metadata/FFQ/FFQ_v10.txt'
ffq = read.table(ffq_file, header = T, sep = '\t')
names(ffq)[1] = 'pn_ID'
ffq = ffq[grepl('A',ffq$pn_ID),]

## 16s
taxa_map_file = '/pita/users/tzipi/projects/multiomics/SOURCE/16S/data/SOURCE_16S_israel_map_v4_noRS.txt'
taxa_map = read.table(taxa_map_file, header = T, sep = '\t')

## mtg
mtg_map_file = '../../metagenomics/data/source_israeli_mtg_species.tsv'
mtg_map = read.table(mtg_map_file, header = T, sep = '\t')

## metadata
metadata_map_file = '/pita/users/tzipi/projects/multiomics/SOURCE/metadata/alona_data/SOURCE_Israel_China_data_v10.txt'
dx_map = read.table(metadata_map_file, header = T, sep = '\t')

pns = unique( c(rna_map$pn_ID, met_map$pn_ID, ffq$pn_ID, taxa_map$pn_ID) )

df = data.frame(pn_ID = pns)
df = merge(df, dx_map[,c('pn_ID','Patient_group2')], by='pn_ID',all.x = T)
# df$Patient_group2[df$pn_ID =='A029'] = 'Israeli_healthy'
df = df[df$pn_ID != 'A029',]
df = df[!is.na(df$Patient_group2),]

df$rnaSeq_TI = df$pn_ID %in% rna_map$pn_ID[rna_map$Location == 'TI'] 
# df$rnaSeq_R = df$pn_ID %in% rna_map$pn_ID[rna_map$Location == 'R'] 

df$metabolomics_stool = df$pn_ID %in% met_map$pn_ID

df$ffq = df$pn_ID %in% ffq$pn_ID

df$taxa_TI1 = df$pn_ID %in% taxa_map$pn_ID[taxa_map$Location_inflamation == 'TI1'] 
df$taxa_TI2 = df$pn_ID %in% taxa_map$pn_ID[taxa_map$Location_inflamation == 'TI2'] 
df$taxa_R4 = df$pn_ID %in% taxa_map$pn_ID[taxa_map$Location_inflamation == 'R4'] 
df$taxa_R3 = df$pn_ID %in% taxa_map$pn_ID[taxa_map$Location_inflamation == 'R3'] 
df$taxa_R3 = df$pn_ID %in% taxa_map$pn_ID[taxa_map$Location_inflamation == 'R3'] 
df$taxa_Stool = df$pn_ID %in% taxa_map$pn_ID[taxa_map$Location_inflamation == 'Stool'] 

df$Metagenomics = df$pn_ID %in% names(mtg_map)

df$Demographics = df$pn_ID %in% dx_map$pn_ID[!is.na(dx_map$Age_years) | !is.na(dx_map$Gender)] 
df$Environment = df$pn_ID %in% dx_map$pn_ID[dx_map$IOIBD_completed == 'Yes'] 
df$CRP = df$pn_ID %in% dx_map$pn_ID[ !is.na(dx_map$CRP_numeric) ]
df$Calprotectin = df$pn_ID %in% dx_map$pn_ID[ !is.na(dx_map$Calprotectin_numeric) ] 


df$fill_sum = rowSums(df[,3:dim(df)[2]])

df$taxa_TI1[df$pn_ID == 'A061'] = F
df$taxa_TI2[df$pn_ID == 'A061'] = T
# df$rnaSeq_TI1[df$pn_ID == 'A061'] = F
# df$rnaSeq_TI2[df$pn_ID == 'A061'] = T

df$Patient_group2[df$Patient_group2 == 'Israeli_Crohns'] = 'Israel CD'
df$Patient_group2[df$Patient_group2 == 'Israeli_healthy'] = 'Israel control'

dfm = reshape::melt(data = df, id.vars = c('pn_ID','fill_sum','Patient_group2'), variable_name = 'Data_type')

dfm$Data_type = gsub('taxa','16S',dfm$Data_type)
dfm$Data_type = gsub('ffq','FFQ',dfm$Data_type)


dfm$pn_ID = factor(dfm$pn_ID, levels = df$pn_ID[order(df$fill_sum)])
dfm$Type = gsub('_.*','',dfm$Data_type)
dfm$Source = gsub('.*_','',dfm$Data_type)

dfm$Type[dfm$Type %in% c('Demographics','Environment','CRP','Calprotectin')] = 'Metadata'
dfm$Type[dfm$Type %in% c('16S','Metagenomics')] = 'Bacteria'
dfm$Type[dfm$Type %in% c('rnaSeq')] = 'RNA-Seq'
dfm$Type[dfm$Type %in% c('metabolomics')] = 'Metabolomics'

dfm$Data_type[dfm$Data_type %in% c('metabolomics_stool')] = 'Metabolomics stool'
dfm$Data_type[dfm$Data_type %in% c('rnaSeq_TI')] = 'RNA-Seq TI'
dfm$Data_type[dfm$Data_type %in% c('rnaSeq_R')] = 'RNA-Seq R'
dfm$Data_type[dfm$Data_type %in% c('16S_TI1')] = '16S TI inflamed'
dfm$Data_type[dfm$Data_type %in% c('16S_TI2')] = '16S TI non-inflamed'
dfm$Data_type[dfm$Data_type %in% c('16S_R3')] = '16S R inflamed'
dfm$Data_type[dfm$Data_type %in% c('16S_R4')] = '16S R non-inflamed'
dfm$Data_type[dfm$Data_type %in% c('Metagenomics')] = 'Metagenomics stool'
dfm$Data_type = gsub('_',' ',dfm$Data_type)


dfm$Type = factor(dfm$Type, levels = c('Metadata','FFQ','Bacteria','RNA-Seq','Metabolomics'))
dfm$Data_type = factor(dfm$Data_type, levels = c('Demographics','Environment','CRP',
                                                 'Calprotectin','FFQ','16S Stool',
                                                 '16S TI inflamed','16S TI non-inflamed',
                                                 '16S R inflamed','16S R non-inflamed',
                                                 'Metagenomics stool','RNA-Seq TI','RNA-Seq R',
                                                 'Metabolomics stool'))

dfm$Patient_group3 = dfm$Patient_group2
for (pg in unique(dfm$Patient_group2))
{
  dfm$Patient_group3[dfm$Patient_group2 == pg] = sprintf('%s\nn=%s', gsub('_',' ',pg), length(unique( dfm$pn_ID[dfm$Patient_group2 == pg] )))
}

# dfm = dfm[dfm$Type !='metabolomics',]

g = ggplot(dfm, aes(x=Data_type, y=pn_ID, alpha = value, fill = Type)) + 
  geom_tile() + 
  # geom_tile(colour = 'gray40') + 
  scale_alpha_manual(values = c(.9,0), name = 'Data\navailable', breaks = c(T,F)) +
  ggtitle('Israel') + 
  scale_fill_manual(values = c('#293462','#006E7F','#F8CB2E','#EE5007','#990000')) + 
  # facet_grid(~Type + Source, scales = 'free') + 
  facet_grid(Patient_group3~. , scales = 'free', space = 'free') + 
  ylab('Subject') + xlab ('')  +
  theme_bw() +scale_y_discrete(expand=c(0,0)) +scale_x_discrete(expand=c(0,0)) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.key = element_rect(colour = "black")) 

# ggsave('Israel_cohort_figure_full.pdf',plot = g, device = 'pdf', width = 4,height = 8, limitsize = FALSE)
ggsave('Israel_cohort_figure_met_full_v2.pdf',plot = g, device = 'pdf', width = 4,height = 8, limitsize = FALSE)

g2 = g+theme(axis.text.y = element_blank(), axis.ticks.y = element_blank() )
# ggsave('Israel_cohort_figure.pdf',plot = g2, device = 'pdf', width = 4,height = 8, limitsize = FALSE)
ggsave('Israel_cohort_figure_met_v2.pdf',plot = g2, device = 'pdf', width = 4,height = 8, limitsize = FALSE)


row.names(dx_map) = dx_map$pn_ID
dx_map = dx_map[df$pn_ID,]
df$Age = dx_map$Age_years
df$Gender = dx_map$Gender
df$BMI = dx_map$BMI
df$CRP_mg_L = dx_map$CRP_numeric
df$Calprotectin_ug_g = dx_map$Calprotectin_numeric

write.table(x = df, file = 'Israel_basic_map.txt', quote = F, sep = '\t',row.names = F)


