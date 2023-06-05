# install.packages(c("matrixStats", "Hmisc", "splines", 
#                    "foreach", "doParallel", "fastcluster", 
#                    "dynamicTreeCut", "survival"))
# # source("http://bioconductor.org/biocLite.R")
# BiocManager::install("GO.db" )
# BiocManager::install("preprocessCore" )
# BiocManager::install("impute" )
# 
# install.packages('/pita/users/tzipi/bin/WGCNA_1.69-81.tar.gz', repos = NULL, lib=.Library)             
library(WGCNA)
options(stringsAsFactors = FALSE);

# source('/pita/users/tzipi/projects/rnaSeq/lncRNA/rnaSeq/R_analysis/lncRNA_funcs.R')

type = 'ProteinCoding' 


# c = c('Protect','SEEM','SOURCE')
c = c('SOURCE')
# country = 'China'
country = 'Israel'
location = 'TI'
# location = 'R'
# location = 'PBMC'
dx = 'all'
# dx = 'Healthy'

oneSample = 'oneSample'
# oneSample = 'allSamples'

corr_control_only_flag = F
if ( corr_control_only_flag )
  corr_disease_instead_of_control_flag = F

name = sprintf('%s_%s_%s_%s_%s', c, country,location, dx, oneSample)
path = 'rnaSeq_main/res/'


if (country == 'China' & location == 'TI' & dx == 'all' & oneSample == 'allSamples')
{
  cut_hight = 50000
  pwr = 9
}

if (country == 'China' & location == 'TI' & dx == 'all' & oneSample == 'oneSample')
{
  cut_hight = 50000
  # cut_hight = 35000
  pwr = 10
}

if (country == 'China' & location == 'TI' & dx == 'Healthy')
{
  cut_hight = 50000
  pwr = 5
}
if (country == 'Israel' & location == 'TI' & dx == 'all')
{
  cut_hight = 70000
  pwr = 10
}
if (country == 'Israel' & location == 'TI' & dx == 'Healthy')
{
  cut_hight = 70000
  pwr = 10
}
if (country == 'Israel' & location == 'PBMC' & dx == 'all')
{
  cut_hight = 70000
  pwr = 10
}

if (country == 'Israel' & location == 'R' & dx == 'all')
{
  cut_hight = 70000
  pwr = 10
}
if (country == 'China' & location == 'R' & dx == 'all')
{
  cut_hight = 40000
  pwr = 8
}

# set needed tables
metadata_file = '../../rnaSeq/China_Israel/data/SOURCE_israel_china_rnaSeq_map_v2.txt'
metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'))
metadata_df$Cohort = 'SOURCE'
metadata_df = metadata_df[metadata_df$Country == country & metadata_df$Location == location,]

## adding Amnon's health index
health_index_file = '../../metadata/health_index_amnon.txt'
health_index = read.table(file = health_index_file, header = T,sep = '\t')
metadata_df$sample_ID = sprintf('%s_%s',metadata_df$pn_ID, metadata_df$location_inflamation)

metadata_df2 = merge(x = metadata_df, y=health_index, by = 'sample_ID', all.x = T, all.y = F)
metadata_df$health_index_same_sample = metadata_df2$health_index

metadata_df3 = merge(x = metadata_df, y=health_index[health_index$location == 'Stool',], by = 'pn_ID', all.x = T, all.y = F)
metadata_df$health_index_stool = metadata_df3$health_index

metadata_df$inflammation_status[metadata_df$Dx == 'healthy'] = 'Non-inflamed'

metadata_df$SampleID =  make.names(metadata_df$SampleID)

TPM_file = sprintf('../../rnaSeq/China_Israel/res/SOURCE_%s_%s_geneFiltered_txi_res.txt', country, location)
tpm = read.table(file = TPM_file, header = T,sep = '\t', stringsAsFactors = F,row.names = 1)

if (dx == 'Healthy') ## filter to healthy samples
{
  metadata_df = metadata_df[metadata_df$Dx == 'healthy',]
}
if ( oneSample == 'oneSample' ) ## filter to one sapmle per patient - if there are 2, keep the inflammed sample.
{
  good_pos = c()
  pns = unique(metadata_df$pn_ID)
  for (pn in pns)
  {
    pos = which(metadata_df$pn_ID == pn)
    if (length(pos) == 1)
      good_pos = c(good_pos, pos)
    if (length(pos) == 2)
    {
      p = which(metadata_df$pn_ID == pn & metadata_df$inflammation_status=='Inflamed')
      if(length(p) == 1)
      {
        good_pos = c(good_pos, p)
      } else
        print(sprintf('more than one inflammed sample for patient %s   problem!', pn))
    }
    if (length(pos) > 2)
      print(sprintf('more than 2 samples for patient %s   problem!', pn))
  }
  metadata_df = metadata_df[good_pos,]
}
metadata_df = metadata_df[order(metadata_df$SampleID),]
tpm = tpm[, make.names(sprintf('%s', metadata_df$SampleID))]

wgcna_res_file = 'rnaSeq_main/res/SOURCE_Israel_TI_all_oneSample_withFFQKcalNorm_gene_module_table.txt'
out_path = 'rnaSeq_main/res/module_corrs'

dir.create(out_path)
wr = read.table(wgcna_res_file, sep = '\t',header = T, comment.char = '', quote="\"")

tpm = as.data.frame(t(tpm))
names(tpm) = make.names(names(tpm))

# type = 'kcal_norm'
# ffq_file = sprintf('../../metadata/FFQ/FFQ_v9_%s.txt',type)
# ffq = read.table(ffq_file, header = T, sep = '\t', row.names = 1)
# 
# row.names(tpm) = gsub('S','A',metadata_df$pn_ID)
# tpm2 = tpm; tpm2$pn_ID = row.names(tpm2)
# ffq$pn_ID = row.names(ffq)
# tpm2 = merge(tpm2, ffq, by='pn_ID',all.x = T)
# row.names(tpm2) = tpm2$pn_ID
# tpm2 = tpm2[gsub('S','A',metadata_df$pn_ID),]
# tpm = tpm2[,names(tpm2)!='pn_ID']

for ( module in unique(wr$module_color) )
{
  wr_f = wr[wr$module_color == module,]
  mm_name = sprintf('MM%s',module)
  wr_f = wr_f[order(wr_f[[mm_name]]),]
  
  cor_type = 'spearman'
  # res <- cor(df2, method = cor_type, na.rm=T)
  library("Hmisc")
  res2 <- rcorr(as.matrix(tpm[,make.names(wr_f$Gene)]), type = cor_type)
  
  col2 = colorRampPalette(rev( c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                                 "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061") ))(200)
  # library(corrplot)
  # corrplot(res2$r, type = "upper", order = "hclust", 
  #          tl.col = "black", method = 'circle', col = col2, tl.cex= 0.8)
  
  library(reshape)
  res_melt = reshape::melt(data = res2$r)
  temp = reshape::melt(data = res2$P)
  res_melt$p_val = temp$value
  res_melt$q_val = p.adjust(res_melt$p_val, method = 'BH', n = length(res_melt$p_val))
  
  data <- t(res2$r)
  ord <- hclust( dist(data, method = "euclidean"), method = "ward.D" )$order
  res_melt$X1 = factor(res_melt$X1, row.names(res2$r)[ord])
  res_melt$X2 = factor(res_melt$X2, row.names(res2$r)[ord])
  
  cor_heatmap = ggplot(res_melt, aes(x=X1, y=X2, fill = value)) + 
    geom_tile() + 
    # scale_fill_gradientn(colours = col2) +
    scale_fill_gradientn(colours = col2, limits = c(-1,1)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    # geom_point(aes(shape = q_val<=0.05)) + scale_shape_manual(values = c(NA, 8)) +
    # geom_text(aes(label=sprintf('%.2f\n%.2e',value,q_val )), size=3)
    scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
    # geom_text(aes(label=sprintf('%.2f',value )), size=3) + 
    xlab('') + ylab('') + ggtitle(module)
  ggsave(sprintf('%s/%s_mets_corrs.jpg', out_path,module),plot = cor_heatmap, device = 'jpg', width = 15,height = 15)
  
}


