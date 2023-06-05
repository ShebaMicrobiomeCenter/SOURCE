source('/pita/users/tzipi/code/R_figs/figs_funcs.R')

library(RColorBrewer)
library(ggplot2)
library(networkD3)

out_path = 'res/sankey/'
dir.create(out_path)

specific_var = 'g__.Ruminococcus..s__gnavus_ASV09357'
specific_name = '_gnavus'
# specific_var = 'g__Oxalobacter.s__formigenes_ASV11831'
# specific_name = '_oxalobacter'
# specific_var = ''
# specific_name = ''

dir = '_pos'
# dir = '_neg'

group = '_rural'

min_q = 0.25

filter_halla_res_to_top = function(halla_res, n = 50, both_sides_flag = F)
{
  if ( !both_sides_flag )
  {
    halla_res = halla_res[order(halla_res$association, decreasing = T),]
    halla_res = halla_res[1:n,]
    return(halla_res)
  } else
  {
    temp = halla_res[halla_res$met_direction == 'module_neg',]
    temp = temp[order(temp$association, decreasing = T),]
    temp = temp[1:min( n, dim(temp)[1]),]
    temp2 = halla_res[halla_res$met_direction == 'module_pos',]
    temp2 = temp2[order(temp2$association, decreasing = T),]
    temp2 = temp2[1:min( n, dim(temp2)[1]),]
    halla_res = rbind(temp, temp2)
  }
  halla_res = halla_res[!is.na(halla_res$X_features),]
  return(halla_res)
}

# ## filter metabolomics data using WGCNA modules
# modules_mets_file = '../WGCNA/rnaSeq_main/res_v2/SOURCE_Israel_TI_withStoolMetablolomicsV2LogClean_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default//selected_modules//SOURCE_Israel_TI_withStoolMetablolomicsV2LogClean_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default_module_cor_df.txt'
# md_mets = read.table(file = modules_mets_file, header = T, sep = '\t')
# md_mets$variable = make.names(md_mets$variable)
# metadata_prms = c('Gender','Age_years','BMI','Disease','inflammed','CRP_numeric','Calprotectin_numeric','Active_Smoker_manual','Previous_or_Current_Tobacco_Use_manual','health_index_same_sample','health_index_stool')
# md_mets = md_mets[! md_mets$variable %in% c(metadata_prms),]
# 
# wanted_modules = c('yellow','green','red','pink', 'purple','tan','salmon','black','brown')
# 
# dx_down_modules = c('yellow','green','red','pink')
# dx_up_modules = c('purple','tan','salmon','black','brown')

col2 = colorRampPalette(rev( c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                               "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061") ))(200)

data_type = 'FFQ'
# data_type = '16S_Stool'
# data_type = 'metagenomics_species'
# data_type = 'metagenomics_pathabundance'
# data_type = 'metagenomics_ecsNamed'

halla_res_all = data.frame()

# data_types = c('metagenomics_species','metagenomics_pathabundance','metagenomics_ecsNamed')
data_types = c('metabolomics_Stool','FFQ','16S_Stool')
data_types2 = c('metabolomics_Stool','FFQ','16S_Stool')
# md = 'tan'
# for ( i in 1:length(wanted_modules) )
for ( data_type in data_types )
{
  for (data_type2 in data_types2 )
  {
    # md = wanted_modules[i]
    
    # md_up_mets = md_mets$variable[ md_mets$q_value <= 0.25 & md_mets$module == md & md_mets$correlation > 0]
    # md_down_mets = md_mets$variable[ md_mets$q_value <= 0.25 & md_mets$module == md & md_mets$correlation < 0]
    
    halla_file = sprintf('res/full_pairs_china_rural/China_%s_vs_China_%s%s/halla_res_alpha0.25/all_associations.txt', data_type, data_type2, group)
    if (file.exists(halla_file))
    {
      halla_res_full = read.table(file = halla_file, header = T, sep = '\t', comment.char = '',quote = '')
    } else
      next()
    # halla_res_full$met_direction = ifelse(halla_res_full$X_features %in% md_up_mets, 'module pos','module neg')
    halla_res = halla_res_full[ halla_res_full$q.values <= 0.25, ]
    if ( dir == '_pos')
      halla_res = halla_res[ halla_res$association >0, ]
    if ( dir == '_neg')
      halla_res = halla_res[ halla_res$association <0, ]
    
    # halla_res = halla_res[halla_res$met_direction == 'module pos',]
    
    if (dim(halla_res)[1]==0)
      next()
    
    links <- data.frame(
      source=halla_res$X_features, 
      target=halla_res$Y_features, 
      value=rep(x = 1, length(halla_res$X_features))
      # value=halla_res$association
    )
    
    # From these flows we need to create a node data frame: it lists every entities involved in the flow
    nodes <- unique( data.frame( name=c(as.character(links$source), 
                                        as.character(links$target)) )
    )
    
    # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
    links$IDsource <- match(links$source, nodes$name)-1 
    links$IDtarget <- match(links$target, nodes$name)-1
    
    p <- sankeyNetwork(Links = links, Nodes = nodes,
                       Source = "IDsource", Target = "IDtarget",
                       Value = "value", NodeID = "name", 
                       sinksRight=FALSE, fontSize = 30,  nodeWidth = 10 )
    p
    
    # library(htmlwidgets)
    # saveWidget(p, file=sprintf('%s/%s_%s_pos.html',out_path, md, data_type))
    # 
    # library(webshot)
    # #install phantom:
    # # webshot::install_phantomjs()
    # # Make a webshot in pdf : high quality but can not choose printed zone
    # webshot("/HtmlWidget/sankeyBasic1.html" , "/HtmlWidget/sankeyBasic1.pdf", delay = 0.2)
    halla_res$data_type = data_type
    halla_res$data_type2 = data_type2
    
    # halla_res = filter_halla_res_to_top(halla_res, n = 25, both_sides_flag =F)
    
    halla_res_all = rbind(halla_res_all, halla_res)
  }
}
halla_res_all = halla_res_all[!is.na(halla_res_all$X_features),]
# halla_res_all$Y_features[halla_res_all$data_type == 'metagenomics_species'] = 
#   gsub('.*[|]','', halla_res_all$Y_features[halla_res_all$data_type == 'metagenomics_species'])

halla_res_all =  halla_res_all[halla_res_all$q.values<=min_q,]
if (specific_var != '')
{
  temp =  halla_res_all[halla_res_all$X_features == specific_var |
                          halla_res_all$Y_features == specific_var ,]
  temp = unique(c(temp$X_features, temp$Y_features))
  halla_res_all = halla_res_all[ halla_res_all$X_features %in% temp & 
                                   halla_res_all$Y_features %in% temp,  ]
}


library(ggalluvial)
g = ggplot(data = halla_res,
           aes(axis1 = X_features, axis2 = Y_features)) +
  # scale_x_discrete(limits = c("Class", "Sex", "Age"), expand = c(.2, .05)) +
  # xlab("Demographic") +
  geom_alluvium() +
  geom_stratum() +
  # facet_grid(data_type~., scale = 'free') + 
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal()  


# create a data frame giving the hierarchical structure
# d1 <- data.frame(from="origin", to=c('Metabolite',unique(halla_res_all$data_type)))
# d2 <- data.frame(from='Metabolite', to=unique(halla_res_all$X_features))
# d3 = data.frame()
# for ( data_type in unique(halla_res_all$data_type) )
# {
#   temp = data.frame(from=data_type, 
#                     to=unique(halla_res_all$Y_features[halla_res_all$data_type == data_type]))
#   d3 = rbind(d3, temp)
# }
d1 <- data.frame(from="origin", to=c(unique(c(halla_res_all$data_type,halla_res_all$data_type2))) )
# d2 <- data.frame(from='Metabolite', to=unique(halla_res_all$X_features))
d2 = data.frame()
for ( data_type2 in unique(halla_res_all$data_type2) )
{
  temp = data.frame(from=data_type2, 
                    to=unique(halla_res_all$Y_features[halla_res_all$data_type2 == data_type2]))
  d2 = rbind(d2, temp)
}
d3 = data.frame()
for ( data_type in unique(halla_res_all$data_type) )
{
  temp = data.frame(from=data_type, 
                    to=unique(halla_res_all$X_features[halla_res_all$data_type == data_type]))
  d3 = rbind(d3, temp)
}
hierarchy <- rbind(d1, d2, d3)
hierarchy = hierarchy[rev(order(hierarchy$from)),]

connect = data.frame(from = halla_res_all$X_features, to = halla_res_all$Y_features, value = halla_res_all$association)

# create a vertices data.frame. One line per object of our hierarchy
vertices  <-  data.frame(
  name = unique(c(as.character(hierarchy$from), as.character(hierarchy$to))) , 
  value = runif(length(unique(c(as.character(hierarchy$from), as.character(hierarchy$to)))))
)
# Let's add a column with the group of each name. It will be useful later to color points
vertices$group  <-  hierarchy$from[ match( vertices$name, hierarchy$to ) ]

#Let's add information concerning the label we are going to add: angle, horizontal adjustement and potential flip
#calculate the ANGLE of the labels
vertices$id <- NA
myleaves <- which(is.na( match(vertices$name, hierarchy$from) ))
nleaves <- length(myleaves)
vertices$id[ myleaves ] <- seq(1:nleaves)
vertices$angle <- 90 - 360 * vertices$id / nleaves

# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
vertices$hjust <- ifelse( vertices$angle < -90, 1, 0)

# flip angle BY to make them readable
vertices$angle <- ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)


# Create a graph object
mygraph <- igraph::graph_from_data_frame( hierarchy, vertices=vertices )

# The connection object must refer to the ids of the leaves:
from  <-  match( connect$from, vertices$name)
to  <-  match( connect$to, vertices$name)
# Basic graph
name2 = vertices$name
temp1 = as.data.frame(table(connect$from))
temp1 = temp1[temp1$Freq==1,]
temp2 = as.data.frame(table(connect$to))
temp2 = temp2[temp2$Freq==1,]
# name2[name2 %in% c(temp1$Var1, temp2$Var1) ] = ''

out_path = 'res/circos/'
library(ggraph)
gg = ggraph(mygraph, layout = 'dendrogram', circular = T) + 
  geom_conn_bundle(data = get_con(from = from, to = to), 
                   alpha=0.7, colour="gray", 
                   # tension =0.2 ) + 
                   tension =0.1 ) + 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour = group), size=3) +
  # scale_colour_manual(values= rep( brewer.pal(9,"Set1") , 30)) + 
  scale_colour_manual(values= c('16S_Stool' = '#E41A1C','FFQ' = '#377EB8','metabolomics_Stool' = '#4DAF4A'), name='') + 
  theme_void() + coord_cartesian(clip = "off") + 
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))

# w=6; h=5
w=4; h=3.3
ggsave(sprintf( '%s/china%s%s%s_circle_q%s_noName.pdf', out_path, group, specific_name, dir, min_q),plot = gg, device = 'pdf', width =w,height = h)

name2 = gsub('_',' ', name2)
# name2 = gsub(' [(]','\n(',name2)
gg_named = gg + geom_node_text(aes(x = x*1.1, y=y*1.1, 
                                   filter = leaf, label=name2, 
                                   angle = angle, hjust=hjust), size=3, alpha=1) +
  theme(
    plot.margin=unit(c(7,5,7,7),"cm"),
    # plot.margin=unit(c(7,3,7,7),"cm"),
    legend.position = 'none'
    # legend.position = c(1.5,-0.5)
  )  
gg_named

# w=7.3; h=8 # ox pos
w=6.3; h=7 # ox neg
w=5.8; h=6.5 # gnavus neg
# w=12; h=12
# w=20; h=20
ggsave(sprintf( '%s/china%s%s%s_circle_q%s.pdf', out_path, group, specific_name, dir, min_q),plot = gg_named, device = 'pdf', width =w,height = h)


