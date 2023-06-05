source('/pita/users/tzipi/code/R_figs/figs_funcs.R')


## read israel china metabolism correlation data
modules_mets_file_israel = 'rnaSeq_main/res_v2/SOURCE_Israel_TI_withStoolMetablolomicsV2LogClean_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default//selected_modules//SOURCE_Israel_TI_withStoolMetablolomicsV2LogClean_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default_module_cor_df.txt'
met_i = read.table(file = modules_mets_file_israel, header = T, sep = '\t')
modules_mets_file_china = 'rnaSeq_main/res_other_modules/TI_china_israel_modules_metsLogClean/selected_modules/SOURCE_China_ageMatch_TI_by_SOURCE_Israel_TI_withFFQv10_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default_module_cor_df.txt'
met_c = read.table(file = modules_mets_file_china, header = T, sep = '\t')

met_i$country = 'Israel'
met_c$country = 'China'
mets_all = rbind(met_i, met_c)

## read table of shared metabolites I got from nina
nina_class_file = '../../metabolomics/China/data/Israel_China_names_overlap_nina.csv'
na_str = c('no_data','_','NA','unknown', 'other','na')
nina_class = read.csv(file = nina_class_file, header = T, na.strings = na_str)

## set data
met_i$variable = make.names(met_i$variable)
met_c$variable = make.names(met_c$variable)

# filter to shared metabolites
met_i = met_i[make.names(met_i$variable) %in% make.names(nina_class$Feature),]
met_c = met_c[make.names(met_c$variable) %in% make.names(nina_class$Feature),]

# filter to significant associations
met_i = met_i[met_i$q_value<=0.25,]
met_c = met_c[met_c$q_value<=0.25,]

# met_i = met_i[met_i$variable %in% met_i$variable[met_i$q_value<=0.25],]
# met_c = met_c[met_c$variable %in% met_c$variable[met_c$q_value<=0.25],]

## merge data
mets = rbind(met_i,met_c)
mets$pair = sprintf('%s_%s',mets$variable, mets$module)

# look at directions
up = unique(mets$pair[mets$q_value<=0.25])
up_res = vector(mode = 'character',length = length(up))
for ( i in 1:length(up) )
{
  p = up[i]
  pos_i = which(mets$pair == p & mets$country == 'Israel')
  pos_c = which(mets$pair == p & mets$country == 'China')
  if ( length(pos_i) == 0 & length(pos_c) == 1 )
  {
    up_res[i] = 'Only china'
  } else if ( length(pos_c) == 0 & length(pos_i) == 1 )
  {
    up_res[i] = 'Only israel'
  } else if ( length(pos_i) == 1 & length(pos_i) == 1 )
  {
    if (mets$correlation[pos_i] > 0 & mets$correlation[pos_c] > 0 )
    {
      up_res[i] = 'both up'
    } else if (mets$correlation[pos_i] < 0 & mets$correlation[pos_c] < 0 )
    {
      up_res[i] = 'both down'
    } else
      up_res[i] = 'different directions'
  } else
    up_res[i] = 'problem'
}

## seperate metabolites to control associated and disease associated, by module disease association
control_mds = c('yellow','green','red','pink')
disease_mds = c('purple','tan','salmon','black','brown')

mets$cor_dir = vector(mode = 'character', length = dim(mets)[1])
mets$cor_dir[mets$module %in% control_mds & mets$correlation > 0 ] = 'Control associated'
mets$cor_dir[mets$module %in% control_mds & mets$correlation < 0 ] = 'Disease associated'
mets$cor_dir[mets$module %in% disease_mds & mets$correlation > 0 ] = 'Disease associated'
mets$cor_dir[mets$module %in% disease_mds & mets$correlation < 0 ] = 'Control associated'

## filter to 1 strongest correaltion per metabolite (clean repeating metabolites)
mets_f = mets
good_pos = c()
for ( mc in unique( sprintf('%s_%s', mets_f$variable, mets_f$country) ) )
for ( i in 1:dim(mets_f)[1] )
{
  pos1 = mets_f$country == mets_f$country[i] & mets_f$variable == mets_f$variable[i]
  max_cor = max( mets_f$correlation[pos1] )
  pos2 = which( mets_f$correlation == max_cor & pos1 )
  if ( length(pos2) == 1 )
  {
    good_pos = c(good_pos, pos2)
  } else
    print ( 'problem!')
    
}
mets_f = mets_f[unique(good_pos),]
mets_f = unique(mets_f[,c('variable','country','cor_dir')])

## set results table
mc = mets_f[mets_f$country == 'China',]
mi = mets_f[mets_f$country == 'Israel',]
mic = merge(mi,mc,by ='variable', all = T)
mic$cor_dir.x[is.na(mic$cor_dir.x)] = 'not sig'
mic$cor_dir.y[is.na(mic$cor_dir.y)] = 'not sig'
mic$both_dir = sprintf('Israel %s China %s', mic$cor_dir.x, mic$cor_dir.y)

table(mic$both_dir)

mic$both_cmp = mic$both_dir
mic$both_cmp[mic$both_cmp %in% c('Israel Control associated China Control associated',
                                 'Israel Disease associated China Disease associated')] = 
  'Same ditrection'
mic$both_cmp[mic$both_cmp %in% c('Israel Control associated China Disease associated',
                                 'Israel Disease associated China Control associated')] = 
  'Different ditrections'
mic$both_cmp[mic$both_cmp %in% c('Israel Control associated China not sig',
                                 'Israel Disease associated China not sig')] = 
  'Only Israel'
mic$both_cmp[mic$both_cmp %in% c('Israel not sig China Control associated',
                                 'Israel not sig China Disease associated')] = 
  'Only China'
table(mic$both_cmp)


# create veanns of results - how many shared matbolites in control and disease asosciated
# install.packages("VennDiagram")
library("VennDiagram")
grid.newpage()                    # Create new plotting page
venn.plot = draw.pairwise.venn(area1 = sum(mic$cor_dir.x == 'Control associated'),    # Draw pairwise venn diagram
                   area2 = sum(mic$cor_dir.y == 'Control associated'),
                   cross.area =  sum(mic$both_dir == 'Israel Control associated China Control associated'),
                   fill = c("light blue", "pink"), 
                   category = c("Israel\nControl", "China \nControl"),
                   # category = c("Israel", "China"),
                   alpha = 0.5, lty = "blank", cat.pos = 180)
out_path = 'rnaSeq_main/israel_china_mets_cmp/'
dir.create(out_path)
# tiff(filename = sprintf('%s/wgcna_israel_china_mets_control_venn.tiff', out_path),
#   compression = "lzw", width = 5, height = 4, units = 'cm', res=300);
# grid.draw(venn.plot);
# dev.off();

grid.newpage()                    # Create new plotting page
venn.plot_d = draw.pairwise.venn(area1 = sum(mic$cor_dir.x == 'Disease associated'),    # Draw pairwise venn diagram
                   area2 = sum(mic$cor_dir.y == 'Disease associated'),
                   cross.area =  sum(mic$both_dir == 'Israel Disease associated China Disease associated'),
                   fill = c("light blue", "pink"), 
                   category = c("Israel\nDisease", "China\nDisease"),
                   alpha = 0.5, lty = rep("blank", 2), cat.pos = 180)
# tiff(filename = sprintf('%s/wgcna_israel_china_mets_disease_venn.tiff', out_path),
#      compression = "lzw", width = 5, height = 4, units = 'cm', res=300);
# grid.draw(venn.plot_d)
# dev.off()

## calculate p value (numbers were calculated manually with Yael and Amnon)
p_israel_ct = 24/39
p_israel_ct = 21/39
prob = p_israel_ct*p_israel_ct + (1-p_israel_ct)*(1-p_israel_ct)
binom.test(x = 32, n = 39 , p = prob)


