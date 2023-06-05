source('/pita/users/tzipi/code/R_figs/figs_funcs.R')

library(RColorBrewer)
library(ggplot2)
library(networkD3)
library(igraph)
library(ggraph)

top_n = 50
filter_flag = T
both_sides_flag = T
extra_name = ifelse(filter_flag,sprintf('top%s',top_n),'')

if ( both_sides_flag & filter_flag ) { extra_name = sprintf('%s_both_sides',extra_name) }

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

out_path = 'res/sankey/'
dir.create(out_path)

## filter metabolomics data using WGCNA modules
modules_mets_file = '../WGCNA/rnaSeq_main/res_v2/SOURCE_Israel_TI_withStoolMetablolomicsV2LogClean_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default//selected_modules//SOURCE_Israel_TI_withStoolMetablolomicsV2LogClean_oneSample_pwr_12_netType_signed_hybrid_minMdSz_30_WGCNA_default_module_cor_df.txt'
md_mets = read.table(file = modules_mets_file, header = T, sep = '\t')
md_mets$variable = make.names(md_mets$variable)
metadata_prms = c('Gender','Age_years','BMI','Disease','inflammed','CRP_numeric',
                  'Calprotectin_numeric','Active_Smoker_manual',
                  'Previous_or_Current_Tobacco_Use_manual','health_index_same_sample','health_index_stool')
md_mets = md_mets[! md_mets$variable %in% c(metadata_prms),]

wanted_modules = c('yellow','green','red','pink', 'purple','tan','salmon','black','brown')

dx_down_modules = c('yellow','green','red','pink')
dx_up_modules = c('purple','tan','salmon','black','brown')

col2 = colorRampPalette(rev( c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                               "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061") ))(200)

data_type = 'FFQ'
# data_type = '16S_Stool'
# data_type = 'metagenomics_species'
# data_type = 'metagenomics_pathabundance'
# data_type = 'metagenomics_ecsNamed'

md = 'black'
for ( i in 1:length(wanted_modules) )
{
  md = wanted_modules[i]
  
  md_up_mets = md_mets$variable[ md_mets$q_value <= 0.25 & md_mets$module == md & md_mets$correlation > 0]
  md_down_mets = md_mets$variable[ md_mets$q_value <= 0.25 & md_mets$module == md & md_mets$correlation < 0]
  
  halla_file = sprintf('res/full_pairs/Israel_metabolomicsM%s_Stool_vs_Israel_%s/halla_res_alpha0.25/all_associations.txt',md, data_type)
  halla_res_full = read.table(file = halla_file, header = T, sep = '\t', comment.char = '',quote = '')
  halla_res_full$met_direction = ifelse(halla_res_full$X_features %in% md_up_mets, 'module_pos','module_neg')
  halla_res = halla_res_full[ halla_res_full$q.values <= 0.25, ]
  halla_res = halla_res[ halla_res$association >0, ]
  # halla_res = halla_res[ halla_res$association <0, ]
  
  if ( filter_flag )
    halla_res = filter_halla_res_to_top(halla_res, n=top_n, both_sides_flag = both_sides_flag)
  
  halla_res$Y_features = gsub('.*[|]','',halla_res$Y_features)
  if (dim(halla_res)[1]==0)
    next()
  
  links <- data.frame(
    source=halla_res$X_features, 
    target=halla_res$Y_features, 
    # value=rep(x = 1, length(halla_res$X_features)),
    value=halla_res$association,
    # group = ifelse(halla_res$met_direction == 'module_pos','1','0')
    group = halla_res$met_direction
  )
  
  # From these flows we need to create a node data frame: it lists every entities involved in the flow
  nodes <- unique( data.frame( name=c(as.character(links$source), 
                                      as.character(links$target)) )
  )
  # nodes$group <- as.factor(c("my_unique_group"))
  # my_color <- 'd3.scaleOrdinal() .domain(["module_neg", "module_pos","my_unique_group"]) .range(["blue", "red","gray"])'
  
  nodes$group <- "target"
  nodes$group[nodes$name %in% links$source[links$group == 'module_pos'] ] = 'pos_cor'
  nodes$group[nodes$name %in% links$source[links$group == 'module_neg'] ] = 'neg_cor'
  if ( md %in% dx_down_modules ) # if those are down modules switch direction, to be disease direction
  {
    my_color <- 'd3.scaleOrdinal() .domain(["module_pos", "module_neg","target","neg_cor","pos_cor"]) .range(["#B4E4FF", "#FFBCBC","#EEEEEE","#F46060","#95BDFF"])'
  } else
    my_color <- 'd3.scaleOrdinal() .domain(["module_neg", "module_pos","target","pos_cor","neg_cor"]) .range(["#B4E4FF", "#FFBCBC","#EEEEEE","#F46060","#95BDFF"])'
  
  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1

  height_num = max(length(unique(links$source)), length(unique(links$target)))
  p <- sankeyNetwork(Links = links, Nodes = nodes,
                     Source = "IDsource", Target = "IDtarget",
                     Value = "value", NodeID = "name", # nodePadding = 5,
                     sinksRight=FALSE, fontSize = 25,  nodeWidth = 10, 
                     # height = 1500, width = 1200
                     # height = 2000, width = 1200
                     colourScale=my_color, LinkGroup="group", NodeGroup="group",
                     height = height_num*26+26, width = 1200
                     # height = 1700, width = 1200
                     # height = 800, width = 1200
                    )
  p
  
  # library(htmlwidgets)
  htmlwidgets::saveWidget(p, file=sprintf('%s/%s_%s_pos_%s.html',out_path, md, data_type, extra_name), )
  # write.table(x = halla_res, file = sprintf('%s/%s_%s_pos.txt',out_path, md, data_type), quote = F, row.names = F, sep = '\t')
   
  # library(webshot)
  # #install phantom:
  # # webshot::install_phantomjs()
  # # Make a webshot in pdf : high quality but can not choose printed zone
  # webshot("/HtmlWidget/sankeyBasic1.html" , "/HtmlWidget/sankeyBasic1.pdf", delay = 0.2)
  
}

# 
# library(ggalluvial)
# g = ggplot(data = halla_res,
#        aes(axis1 = X_features, axis2 = Y_features)) +
#   # scale_x_discrete(limits = c("Class", "Sex", "Age"), expand = c(.2, .05)) +
#   # xlab("Demographic") +
#   geom_alluvium() +
#   geom_stratum() +
#   geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
#   theme_minimal() 
# 
# 
# # create a data frame giving the hierarchical structure
# d1 <- data.frame(from="origin", to=c('Metabolite',data_type))
# d2 <- data.frame(from='Metabolite', to=unique(halla_res$X_features))
# d3 <- data.frame(from=data_type, to=unique(halla_res$Y_features))
# hierarchy <- rbind(d1, d2, d3)
# 
# connect = data.frame(from = halla_res$X_features, to = halla_res$Y_features, value = halla_res$association)
# 
# # create a vertices data.frame. One line per object of our hierarchy
# vertices  <-  data.frame(
#   name = unique(c(as.character(hierarchy$from), as.character(hierarchy$to))) , 
#   value = runif(length(unique(c(as.character(hierarchy$from), as.character(hierarchy$to)))))
# )
# # Let's add a column with the group of each name. It will be useful later to color points
# vertices$group  <-  hierarchy$from[ match( vertices$name, hierarchy$to ) ]
# 
# #Let's add information concerning the label we are going to add: angle, horizontal adjustement and potential flip
# #calculate the ANGLE of the labels
# vertices$id <- NA
# myleaves <- which(is.na( match(vertices$name, hierarchy$from) ))
# nleaves <- length(myleaves)
# vertices$id[ myleaves ] <- seq(1:nleaves)
# vertices$angle <- 90 - 360 * vertices$id / nleaves
# 
# # calculate the alignment of labels: right or left
# # If I am on the left part of the plot, my labels have currently an angle < -90
# vertices$hjust <- ifelse( vertices$angle < -90, 1, 0)
# 
# # flip angle BY to make them readable
# vertices$angle <- ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)
# 
# 
# # Create a graph object
# mygraph <- graph_from_data_frame( hierarchy, vertices=vertices )
# 
# # The connection object must refer to the ids of the leaves:
# from  <-  match( connect$from, vertices$name)
# to  <-  match( connect$to, vertices$name)
# # Basic graph
# gg = ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
#   geom_conn_bundle(data = get_con(from = from, to = to), 
#                    alpha=0.5, colour="skyblue", tension =0.2) + 
#   geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour = group), size=3) +
#   scale_colour_manual(values= rep( brewer.pal(9,"Paired") , 30)) + 
#   theme_void() + 
#   geom_node_text(aes(x = x*1.1, y=y*1.1, 
#                      filter = leaf, label=name, 
#                      angle = angle, hjust=hjust), size=1.5, alpha=1) +
#   theme(
#     plot.margin=unit(c(0,0,0,0),"cm"),
#   ) +
#   expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))
# 
