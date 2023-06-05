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

## uses WGCNA goodSamplesGenes function, printf the names of the removed samples and features in there are any
## returns a filtered ftr_df (=input feature data frame)
WGCNA_initial_QC = function(ftr_df)
{
  gsg = WGCNA::goodSamplesGenes(ftr_df, verbose = 3);
  gsg$allOK
  if (!gsg$allOK)
  {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0) 
      printFlush(paste("Removing genes:", paste(names(ftr_df)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0) 
      printFlush(paste("Removing samples:", paste(rownames(ftr_df)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    ftr_df = ftr_df[gsg$goodSamples, gsg$goodGenes]
  }
  return(ftr_df)
}

## clusters a sample tree and cuts it at a predefined cutoff. used to remove outlier samples
WGCNA_sample_tree = function(ftr_df, cut_hight, path, name)
{
  sampleTree = hclust(dist(ftr_df), method = "average");
  # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
  # The user should change the dimensions if the window is too large or too small.
  sizeGrWindow(12,9)
  #pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  pdf(sprintf('%s/%s_sample_dendro_bofore_cut.pdf',path, name),width = 13,height = 12 )
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
       cex.axis = 1.5, cex.main = 2)
  
  # Plot a line to show the cut
  abline(h = cut_hight, col = "red");
  dev.off()
  # Determine cluster under the line
  clust = WGCNA::cutreeStatic(sampleTree, cutHeight = cut_hight, minSize = 10)
  table(clust)
  # clust 1 contains the samples we want to keep.
  keepSamples = (clust==1)
  
  WGCNA::collectGarbage();
  return(keepSamples)
}

## re-clusters a sample tree with clean sample list, and plot with metadata
WGCNA_sample_tree_with_metadata = function(ftr_df, metadata_df, path, name)
{
  # Re-cluster samples
  sampleTree2 = hclust(dist(ftr_df), method = "average")
  # Convert traits to a color representation: white means low, red means high, grey means missing entry
  traitColors = WGCNA::numbers2colors(metadata_df, signed = FALSE);
  # Plot the sample dendrogram and the colors underneath.
  
  pdf(sprintf('%s/%s_sample_dendro.pdf',path, name),width = 13,height = 10 + dim(metadata_df)[2]/5 )
  WGCNA::plotDendroAndColors(sampleTree2, traitColors,
                      groupLabels = names(metadata_df), 
                      main = "Sample dendrogram and trait heatmap")
  dev.off()
}

## create a plot for power selection
WGCNA_power_options_plot = function(ftr_df)
{
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  # Call the network topology analysis function
  sft = WGCNA::pickSoftThreshold(ftr_df, powerVector = powers, verbose = 5)
  # Plot the results:
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  return(sft)
}

WGCNA_net = function(ftr_df, sft, automatic_power_flag, pwr = NA, network_type = 'signed hybrid', minModuleSize = 30)
{
  # sft$fitIndices
  if ( automatic_power_flag )
  {
    print(sprintf('used estimated power: %s', sft$powerEstimate))
    pwr = sft$powerEstimate
  } else
    print(sprintf('estimated power: %s, used power: %s', sft$powerEstimate, pwr))
  
  net = WGCNA::blockwiseModules(ftr_df, power = pwr, maxBlockSize = 20000,# TOMType = "signed", 
                         networkType = network_type,
                         minModuleSize = minModuleSize,# minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = F,
                         # saveTOMFileBase = "femaleMouseTOM", 
                         verbose = 3)
  
  net_name = sprintf('pwr_%s_netType_%s_minMdSz_%s', pwr, network_type, minModuleSize)
  return(list( net,net_name ) )
}

## plot a deprogram of the net
WGCNA_plot_net_dendro = function(net, path, name)
{
  # open a graphics window
  sizeGrWindow(12, 9)
  # Convert labels to colors for plotting
  mergedColors = labels2colors(net$colors)
  # Plot the dendrogram and the module colors underneath
  pdf(sprintf('%s/%s_cluster_dendrogram.pdf',path, name),width = 10,height = 5 )
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
}

WGCNA_plot_genes_correlaitons_by_module = function(ftr_df, moduleColors, path,cor_type = 'spearman')
{
  path2 = sprintf('%s/genes_correlation_by_module/', path)
  dir.create(path2)
  col2 = colorRampPalette(rev( c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                                 "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061") ))(200)
  for ( module in unique(moduleColors) )
  {
    ftr_df_f = ftr_df[,moduleColors == module]
    
    library("Hmisc")
    res2 <- rcorr(as.matrix(ftr_df_f), type = cor_type)
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
      xlab('') + ylab('') + ggtitle(sprintf('%s (n=%s)',module,sum(moduleColors == module)))
    
    if ( sum(moduleColors == module) > 50 )
    {
      cor_heatmap = cor_heatmap + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                                        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
    }
    ggsave(sprintf('%s/%s_corrs.jpg', path2,module),plot = cor_heatmap, device = 'jpg', width = 6,height = 6)
    
  }
}

clean_metadata_eampty_noVar = function(metadata_df, verbose = T)
{
  ## remove metadata columns with no data or no variance
  for (i in dim(metadata_df)[2]:1)
    if ( sum(!is.na(metadata_df[,i])) <= 1 | var(metadata_df[,i], na.rm = T) ==0 ) 
    { 
      if ( verbose )
        print(sprintf('remove %s from metadata',names(metadata_df)[i]))
      metadata_df = metadata_df[,-i] 
    }
  return(metadata_df)
}

WGCNA_calculate_module_eigengenes_and_corrs = function(ftr_df,moduleColors, metadata_df, cor_type = 'WGCNA_default', name)
{
  MEs = WGCNA::orderMEs( WGCNA::moduleEigengenes(ftr_df, moduleColors)$eigengenes )
  if ( cor_type == 'WGCNA_default')
  { ## use the default pearson + Student p value
    moduleTraitCor = WGCNA::cor(MEs, metadata_df, use = "p");
    moduleTraitPvalue = WGCNA::corPvalueStudent(moduleTraitCor, nrow(ftr_df));
  } else if (cor_type == 'spearman')
  { ## calculate spearman correlation instead
    # calculate default just to use as a data frame template
    moduleTraitCor = WGCNA::cor(MEs, metadata_df, use = "p");
    moduleTraitPvalue =  WGCNA::corPvalueStudent(moduleTraitCor, nrow(ftr_df));
    for (i in 1:dim(MEs)[2])
    {
      for ( j in 1:dim(metadata_df)[2]) 
      {
        res = cor.test(x=MEs[,i], y= metadata_df[,j],method = 'spearman')
        moduleTraitCor[i,j] = res$estimate
        moduleTraitPvalue[i,j] = res$p.value
      }
    }
    # name = sprintf('%s_spearman', name)
  }
  return(list( MEs, moduleTraitCor, moduleTraitPvalue ))
}

WGCNA_ME_met_heatmap_original = function(moduleTraitCor, moduleTraitPvalue, name, path)
{
  ## plot correlation to all samples
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3.1, 3));
  # Display the correlation values within a heatmap plot
  wd = 8
  if (grepl('FFQ', name))
    wd = 30
  pdf(sprintf('%s/%s_module_trait.pdf',path, name),width = wd,height = 6 )
  WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = gsub(pattern = '_per_1000kcal_day',replacement = '_nrm',x = colnames(moduleTraitCor)),
                 yLabels = rownames(moduleTraitCor),
                 ySymbols =  rownames(moduleTraitCor),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 xLabelsAdj = 1,
                 cex.lab=0.9,
                 xLabelsAngle = 30,
                 main = sprintf('%s_module_trait',name) )
  dev.off()
}

WGCNA_ME_met_heatmap_original_prm_related = function(moduleTraitCor, moduleTraitPvalue, name, prm, cutoff, path)
{
  good_cols = row.names(moduleTraitPvalue)[ moduleTraitPvalue[,prm] <= cutoff ]
  if ( length(good_cols >1))
    WGCNA_ME_met_heatmap_original(moduleTraitCor[good_cols,], moduleTraitPvalue[good_cols,], sprintf('%s_%sF%s',name,prm, cutoff), path)
}

WGCNA_ME_met_heatmap = function(moduleTraitCor, moduleTraitPvalue, name, metadata_prms = c(), min_p = 1, wanted_modules_order = c(), fdr_flag_non_metadata = F, cluster_cols_order_flag = T, show_p_value_fdr_flag = F)
{
  res = WGCNA_ME_met_heatmap_full_res(moduleTraitCor, moduleTraitPvalue, name, metadata_prms, min_p, wanted_modules_order, fdr_flag_non_metadata, cluster_cols_order_flag, show_p_value_fdr_flag)
  return(res[[1]])
}

WGCNA_ME_met_heatmap_full_res = function(moduleTraitCor, moduleTraitPvalue, name, metadata_prms = c(), min_p = 1, wanted_modules_order = c(), fdr_flag_non_metadata = F, cluster_cols_order_flag = T, show_p_value_fdr_flag = F)
{
  if ( !fdr_flag_non_metadata )
  {
    good_pos = apply(moduleTraitPvalue,2,min) <= min_p
    good_pos[colnames(moduleTraitPvalue) %in% metadata_prms] = T
    
    moduleTraitCor = moduleTraitCor[, good_pos]
    moduleTraitPvalue = moduleTraitPvalue[, good_pos]
  }
  if ( min_p != 1)
    name = sprintf('%s\nP <= %s in at least 1 modules',name, min_p)
  
  mdc = reshape::melt(moduleTraitCor)
  temp = reshape::melt(moduleTraitPvalue)
  names(mdc) = c('module','variable','correlation')
  mdc$p_value = temp$value
  mdc$module = gsub('ME','',mdc$module)
  
  mdc_all = mdc
  
  if (fdr_flag_non_metadata) ## add fdr on non-metadata varibale (not in metadata_prms)
  {
    pos = ! mdc$variable %in%  metadata_prms
    mdc$q_value = mdc$p_value
    mdc$q_value[pos] = p.adjust(mdc$q_value[pos], method = 'BH', n = length(mdc$q_value[pos]))
    
    mdc_all = mdc
    ## re-filter by q. value cutoff
    good_vrs = c()
    for ( v in unique(mdc$variable)  )
    {
      if (  min( mdc$q_value[ mdc$variable == v ] ) <= min_p )
      {
        good_vrs = c(good_vrs, v)
      }
    }
    good_vrs = c(good_vrs, metadata_prms[metadata_prms %in% mdc$variable] )
    ## filter data frames to significant variables only
    mdc = mdc[mdc$variable %in% good_vrs,]
    moduleTraitCor = moduleTraitCor[,good_vrs]
    moduleTraitPvalue = moduleTraitPvalue[,good_vrs]
  }
  
  col2 = colorRampPalette(rev( c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                                 "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061") ))(200)
  
  
  temp = moduleTraitCor[row.names(moduleTraitCor) !='MEgrey',]
  if ( length(wanted_modules_order) == 0)
  {
    ord <- hclust( dist(temp, method = "euclidean"), method = "ward.D" )$order
    mdc$module = factor(mdc$module, c(gsub('ME','',row.names(temp))[ord],'grey'))
  } else 
  {
    wanted_modules_order = gsub('ME','',wanted_modules_order)
    if ( ! all(unique(mdc$module) %in% wanted_modules_order))
      print('wanted_modules_order does not include all modules')
    else
      mdc$module = factor(mdc$module, wanted_modules_order)
  }
  good_cols = colnames(moduleTraitCor)[!colnames(moduleTraitCor) %in% metadata_prms]
  temp = moduleTraitCor[,good_cols]
  if ( length(good_cols) > 1 ) # if we have enough data to cluster on
  {
    if ( cluster_cols_order_flag == T )
    {
      ord <- hclust( dist(t(temp), method = "euclidean"), method = "ward.D" )$order
      mdc$variable = factor(mdc$variable, c(colnames(temp)[ord],rev(metadata_prms)) )
    } else
    {
      mdc$variable = factor(mdc$variable, c(colnames(temp),rev(metadata_prms)) )
    }
    
  } else 
  {
    mdc$variable = factor(mdc$variable, c(good_cols,rev(metadata_prms)) )
  }
  
  
  library(ggplot2)
  heatmap = ggplot(mdc, aes(x=module,y=variable,fill = correlation)) +
    geom_tile() + 
    scale_fill_gradientn(colours = col2, limits = c(-1,1)) +
    # geom_point(aes(shape = p_value<=0.05)) + scale_shape_manual(values = c(NA, 8)) +
    # geom_text(aes(label=sprintf('%.2f\n%.2e',correlation,p_value )), size=3) + 
    scale_colour_manual(values = c('black','white')) + 
    # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size= 10)) + 
    theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
    # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size= 10)) + coord_flip() +
    scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
    ylab('') + ggtitle(name)
  if ( (show_p_value_fdr_flag  & fdr_flag_non_metadata) | !fdr_flag_non_metadata )
  {
    heatmap = heatmap + geom_text(aes(label=signif(p_value, 1), colour =p_value < 0.001 ), size=3) 
  } else
    heatmap = heatmap + geom_text(aes(label=signif(q_value, 1), colour =p_value < 0.001 ), size=3) 
  return( list(heatmap, mdc,mdc_all) )
}

WGCNA_ME_met_heatmap_pq = function(moduleTraitCor, moduleTraitPvalue, name, metadata_prms = c(), 
                                   min_p = 1, wanted_modules_order = c(), # fdr_flag_non_metadata = F, 
                                   cluster_cols_order_flag = T, show_p_value_fdr_flag = T, min_q = 1)
{
  res_p = WGCNA_ME_met_heatmap_full_res(moduleTraitCor[wanted_modules,], moduleTraitPvalue[wanted_modules,], name, 
                                        metadata_prms, min_p = min_p, wanted_modules_order = wanted_modules, 
                                        fdr_flag_non_metadata = F, show_p_value_fdr_flag = show_p_value_fdr_flag )
  res_q = WGCNA_ME_met_heatmap_full_res(moduleTraitCor[wanted_modules,], moduleTraitPvalue[wanted_modules,], name, 
                                        metadata_prms, min_p = min_q, wanted_modules_order = wanted_modules, 
                                        fdr_flag_non_metadata = T, show_p_value_fdr_flag = show_p_value_fdr_flag )
  mdc_p = res_p[[2]]
  mdc_q = res_q[[3]]
  mdc_q = mdc_q[mdc_q$variable %in% mdc_p$variable,]
  mdc = mdc_q
  
  temp = moduleTraitCor[row.names(moduleTraitCor) !='MEgrey',]
  if ( length(wanted_modules_order) == 0)
  {
    ord <- hclust( dist(temp, method = "euclidean"), method = "ward.D" )$order
    mdc$module = factor(mdc$module, c(gsub('ME','',row.names(temp))[ord],'grey'))
  } else 
  {
    wanted_modules_order = gsub('ME','',wanted_modules_order)
    if ( ! all(unique(mdc$module) %in% wanted_modules_order))
      print('wanted_modules_order does not include all modules')
    else
      mdc$module = factor(mdc$module, wanted_modules_order)
  }
  good_cols = colnames(moduleTraitCor)[!colnames(moduleTraitCor) %in% metadata_prms & colnames(moduleTraitCor) %in% mdc$variable]
  temp = moduleTraitCor[,good_cols]
  if ( length(good_cols) > 1 ) # if we have enough data to cluster on
  {
    if ( cluster_cols_order_flag == T )
    {
      ord <- hclust( dist(t(temp), method = "euclidean"), method = "ward.D" )$order
      mdc$variable = factor(mdc$variable, c(colnames(temp)[ord],rev(metadata_prms)) )
    } else
    {
      mdc$variable = factor(mdc$variable, c(colnames(temp),rev(metadata_prms)) )
    }
    
  } else 
  {
    mdc$variable = factor(mdc$variable, c(good_cols,rev(metadata_prms)) )
  }
  library(ggplot2)
  col2 = colorRampPalette(rev( c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                                 "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061") ))(200)
  heatmap = ggplot(mdc, aes(x=module,y=variable,fill = correlation)) +
    geom_tile() + 
    geom_tile(data = mdc[mdc$q_value <= min_q & !mdc$variable %in% metadata_prms,], 
              color = 'black', size = 1, aes(x=module,y=variable,fill = correlation)) +
    scale_fill_gradientn(colours = col2, limits = c(-1,1)) +
    # geom_point(aes(shape = p_value<=0.05)) + scale_shape_manual(values = c(NA, 8)) +
    # geom_text(aes(label=sprintf('%.2f\n%.2e',correlation,p_value )), size=3) + 
    scale_colour_manual(values = c('black','white')) + 
    # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size= 10)) + 
    theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
    # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size= 10)) + coord_flip() +
    scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
    ylab('') + ggtitle(name)
  # if ( (show_p_value_fdr_flag  & fdr_flag_non_metadata) | !fdr_flag_non_metadata )
  # {
  heatmap = heatmap + geom_text(aes(label=signif(p_value, 1), colour =p_value < 0.001 ), size=3) 
  # } else
  # heatmap = heatmap + geom_text(aes(label=signif(q_value, 1), colour =p_value < 0.001 ), size=3) 
  return(list(heatmap, mdc))
}

WGCNA_ME_met_heatmap_prm_related = function(moduleTraitCor, moduleTraitPvalue, name, prm, cutoff, metadata_prms = c())
{
  good_cols = row.names(moduleTraitPvalue)[ moduleTraitPvalue[,prm] <= cutoff ]
  if ( length(good_cols >1))
    WGCNA_ME_met_heatmap(moduleTraitCor[good_cols,], moduleTraitPvalue[good_cols,], sprintf('%s_%sF%s',name,prm, cutoff), metadata_prms)
}


WGCNA_ME_met_heatmap_prm_related_pflt = function(moduleTraitCor, moduleTraitPvalue, name, prm, cutoff = 0.05, metadata_prms = c(), p_cutoff = 0.05)
{
  good_cols = row.names(moduleTraitPvalue)[ moduleTraitPvalue[,prm] <= cutoff ]
  good_rows = which( apply(moduleTraitPvalue[good_cols,],2,min) <= p_cutoff | colnames(moduleTraitPvalue) %in% metadata_prms )
  if ( length(good_cols >1) )
    WGCNA_ME_met_heatmap(moduleTraitCor[good_cols,good_rows], moduleTraitPvalue[good_cols,good_rows], 
                         sprintf('%s_%sF%s_pF%s',name,prm, cutoff, p_cutoff), metadata_prms)
}

WGCNA_set_results_table = function(ftr_df, metadata_df, MEs, disease_prm = 'Disease')
{
  # Define variable Disease containing the weight column of datTrait
  Disease = as.data.frame(metadata_df[[disease_prm]]);
  names(Disease) = "Disease"
  # names (colors) of the modules
  modNames = substring(names(MEs), 3)
  
  geneModuleMembership = as.data.frame(WGCNA::cor(ftr_df, MEs, use = "p"));
  MMPvalue = as.data.frame(WGCNA::corPvalueStudent(as.matrix(geneModuleMembership), nrow(ftr_df)));
  
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  
  geneTraitSignificance = as.data.frame(cor(ftr_df, Disease, use = "p"));
  GSPvalue = as.data.frame(WGCNA::corPvalueStudent(as.matrix(geneTraitSignificance), nrow(ftr_df)));
  
  names(geneTraitSignificance) = paste("GS.", names(Disease), sep="");
  names(GSPvalue) = paste("p.GS.", names(Disease), sep="");
  
  ADJ1=abs(cor(ftr_df,use="p"))^pwr
  Alldegrees=WGCNA::intramodularConnectivity(ADJ1, moduleColors)
  
  df = data.frame(row.names = names(ftr_df), module_color = moduleColors)
  
  df = merge(x = df, y = Alldegrees, by = 0); rownames(df)=df$Row.names; df=df[2:dim(df)[2]];
  df = merge(x = df, y = geneTraitSignificance, by = 0); rownames(df)=df$Row.names; df=df[2:dim(df)[2]];
  df = merge(x = df, y = GSPvalue, by = 0); rownames(df)=df$Row.names; df=df[2:dim(df)[2]];
  df = merge(x = df, y = geneModuleMembership, by = 0); rownames(df)=df$Row.names; df=df[2:dim(df)[2]];
  df = merge(x = df, y = MMPvalue, by = 0); rownames(df)=df$Row.names; df=df[2:dim(df)[2]];
  
  df$Gene = row.names(df)
  l = dim(df)[2]
  df = df[,c(l,1:(l-1))]
  
  return(df)
}
