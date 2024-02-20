library(reshape)
library(ggplot2)
# library(easyGgplot2)
library(stringr)
library(matrixStats)
library(vegan)
source('funcs/String_class.R')

# writes only map values(no taxa) into a text file
write_map = function(data_frame, map_name)
{
  ids = names(data_frame);
  no_taxa_ids_pos = !grepl('__',ids) 
  data_frame = data_frame[,no_taxa_ids_pos]
  write.table(data_frame,map_name,sep = "\t" ,quote = F,row.names = F)
}

# writes the full taxa dataframe into a text file
write_taxa = function(data_frame, file_name, transpose_flag = F)
{
  if ( !transpose_flag )
    write.table(data_frame,file_name,sep = "\t" ,quote = F,row.names = F)
  else # write transposed
  { 
    write.table(t(data_frame),file_name,sep = "\t" ,quote = F,col.names = F)
  }
  
}


# setting data with metadata for pcoa ploting.
pcoa_set_data = function(pcoa_file, metadata_file, metadata_cols, cols_num = 3, met_na_str = c('no_data','_','NA','unknown', 'other'))
{
  ## reading pcoa
  # read the pcoa values
  pcoa = read.table(pcoa_file,sep="\t", header=F,skip = 9,fill = T,blank.lines.skip = T)
  pcoa = pcoa[ 1:( dim(pcoa)[1]-2), ] # removing last unrelevant lines
  pcoa = pcoa[,1:(1+cols_num)]
  
  # read dimentions
  con<-file(pcoa_file)
  open(con)
  dims = read.table(con,skip=4,nrow=1) #5-th line
  close(con)
  dims = dims[1:cols_num]
  
  ## reading taxa for metadata
  taxa = read.table(metadata_file,sep="\t", header=TRUE, na.strings = met_na_str, comment.char = '', quote="\"")

  ## adding metadata to pcoa
  met = taxa[metadata_cols]
  
  pcoa = data_add_metadata(pcoa, met, 'V1')
  
  return( list(pcoa, dims) )
}

# setting colors by wanted variable for ploting. used for pcoa 3D scatterplot
set_color = function(col_data, cols, data_wanted_levels = NULL)
{
  # if I want more level (colors) than in the data, mostly to fit colors of another figure
  if ( is.null(data_wanted_levels) )
    data_lvs = levels(col_data)
  else
    data_lvs = data_wanted_levels
  # check that the number of colors fit
  if ( length(data_lvs) != length(cols) )
    print('Color number does not fit.')
  else
  {
    temp = rep(cols[1], length(col_data))
    for ( i in 1:length(data_lvs) )
    {
      temp[ grepl( as.character(data_lvs[i]), as.character(col_data) ) ] = cols[i]
    }
    return(temp)
  }
}

# removes taxa with low aboundece ( = less than high_num have a higher frequncy than cutoff) so that the colours will look better
rmv_low_taxa = function(data, taxa_level = 2, cutoff = 0.01, high_num = 10)
{
  # removing unwanted taxa levels
  data = filter_data(data,list(),level = taxa_level )
  all_ids = names(data)
  taxa_ids = all_ids[get_taxa_level_pos( data, taxa_level )]
  new_data = data.frame(SampleID = data$SampleID)
  for ( id in all_ids )
  {
    if ( id %in% taxa_ids)
    {
      if ( length( data[[id]][ data[[id]] > cutoff ] ) > high_num )
        new_data[[id]] = data[[id]]
    } else
      new_data[[id]] = data[[id]]
  }
  return(new_data)
}

# changing low abundance taxa to "other"
rmv_low_taxa_v2 = function(data, taxa_level = 2, high_taxa_num = 6)
{
  # removing unwanted taxa levels
  data = filter_data(data,list(),level = taxa_level )
  all_ids = names(data)
  taxa_ids = all_ids[get_taxa_level_pos( data, taxa_level )]
  taxa_mean_abd = as.numeric(colMeans( data[,taxa_ids ]) )
  high_taxa = taxa_ids[ order(taxa_mean_abd, decreasing = T)][1:high_taxa_num ]
  low_taxa = taxa_ids[ !taxa_ids %in% high_taxa ]
  other = as.numeric(rowSums( data[,low_taxa ] ) )
  data$k__other = other
  good_cols = c( all_ids[ !all_ids %in% taxa_ids], high_taxa, 'k__other'  )
  
  new_data = data[,good_cols]
  return(new_data)
}

# removes taxa with low aboundece ( = less than high_num have a higher frequncy than cutoff) and sums them into "other"
low_taxa_to_other = function(data, taxa_level = 2, cutoff = 0.01, high_num = 10)
{
  # removing unwanted taxa levels
  data = filter_data(data,list(),level = taxa_level )
  all_ids = names(data)
  taxa_ids = all_ids[get_taxa_level_pos( data, taxa_level )]
  new_data = data.frame(SampleID = data$SampleID)
  other = vector(mode = 'numeric',length = dim(data)[1])
  for ( id in all_ids )
  {
    if ( id %in% taxa_ids)
    {
      if ( length( data[[id]][ data[[id]] > cutoff ] ) > high_num )
      {
        new_data[[id]] = data[[id]]
      } else
      {
        other = other + data[[id]]
      }
    } else
      new_data[[id]] = data[[id]]
  }
  return(new_data)
}

# removes functions (from picrust) with low aboundece ( = less than high_num have a higher frequncy than cutoff) so that the colours will look better
rmv_low_function = function(data, cutoff = 0.01, high_num = 10)
{
  all_ids = names(data)
  # get only position the are nor metadata
  pos = which(all_ids == 'Description')
  taxa_ids = all_ids[(pos+1):length(all_ids)]
  new_data = data.frame(SampleID = data$SampleID)
  for ( id in all_ids )
  {
    if ( id %in% taxa_ids)
    {
      if ( length( data[[id]][ data[[id]] > cutoff ] ) > high_num )
        new_data[[id]] = data[[id]]
    } else
      new_data[[id]] = data[[id]]
  }
  return(new_data)
}

# adds metadta to data by shared id column. for data id=X, for metadata id = SampleID.
data_add_metadata = function(data, metadata, ids_name = 'X')
{ 
  # sorting to same order
  metadata= metadata[order(metadata$SampleID),]
  data = data[order(data[[ids_name]]),]
  
  # if there is more metadata for which there is no data
  if (dim(metadata)[1] > dim(data)[1])
  {
    metadata = metadata[metadata$SampleID %in% data[[ids_name]], ]
  }
  
  cn = colnames(metadata)
  # 1 is SampleID which is same in data and was used to compare order
  for ( i in 2:length(cn) )
  {
    data[cn[i]] = metadata[cn[i]]
  }
  return(data)
}

# returns a data frame without specific value in wanted parameter. removes NA in said parameter as well.
na_rmv_by_col = function(data, parm, rm_prm)
{
  data2 = data

  data2 = data2[ !data2[parm] == rm_prm & !is.na( data2[parm] ), ]

  return(data2)
}


# returns a data frame containing only samples with wanted parameters from data by parms 
# gets a data frame (data) and a list (parms) with names that are same as in data and parameters to filter by
filter_data = function(data, parms, levell = NA, melt_flag = F, rmv_na = TRUE)
{
  data2 = data
  cols = names(parms)
  
  for ( i in 1:length(cols) )
  {
    parm = unlist( parms[cols[i]] )
    if ( length(parm) == 2 ) # in range, >= first and < second.
    {
      if (rmv_na)
        data2 = data2[ data2[cols[i]] >= parm[1] & data2[cols[i]] < parm[2] & !is.na(data2[cols[i]]), ]
      else
        data2 = data2[ data2[cols[i]] >= parm[1] & data2[cols[i]] < parm[2], ]
    }
    if ( length(parm) == 1 ) # equals
    {
      if ( !is.na(parm)  ) # can use NA for variables I dont want to check
      {
        if (rmv_na)
          data2 = data2[ data2[cols[i]] == parm & !is.na(data2[cols[i]]), ]
        else
          data2 = data2[ data2[cols[i]] == parm, ]
      }
    }
  }
  
  if ( !is.na(levell) )
  {
    level_taxa_pos = get_taxa_level_pos( data2, levell )
    ids = names(data2);
    no_taxa_ids_pos = !grepl('__',ids) 
    # get only wanted level
    data2 = data2[,no_taxa_ids_pos | level_taxa_pos]
  }
  if (melt_flag)
  {
    ids = names(data2);
    data2 = melt(data2, id = ids[ !grepl('__',ids) ], variable_name = 'taxa', na.rm = F)
  }
  data2 = droplevels(data2)
  return(data2)
}


# returns a data frame containing only samples without wanted parameters from data by parms 
# gets a data frame (data) and a list (parms) with names that are same as in data and parameters to filter by
filter_data_rm = function(data, parms, level = NA, melt_flag = F, rmv_na = TRUE)
{
  if ( length(parms) > 0)
  {
    data2 = data
    cols = names(parms)
    for ( i in 1:length(cols) )
    {
      parm = unlist( parms[cols[i]] )
      for ( j in 1:length(parm) )
      {
        
        if ( rmv_na )
          data2 = data2[ data2[cols[i]] != parm[j] & !is.na(data2[cols[i]]), ]
        else
          data2 = data2[ data2[cols[i]] != parm[j], ]
      }
    }
    
    if ( !is.na(level) )
    {
      level_taxa_pos = get_taxa_level_pos( data2, level )
      ids = names(data2);
      no_taxa_ids_pos = !grepl('__',ids) 
      # get only wanted level
      data2 = data2[,no_taxa_ids_pos | level_taxa_pos]
    }
    if (melt_flag)
    {
      ids = names(data2);
      data2 = melt(data2, id = ids[ !grepl('__',ids) ], variable_name = 'taxa', na.rm = F)
    }
    data2 = droplevels(data2)
  }
  return(data2)
}

# returns a data frame without bad_val (or NA) in  filter_parm. can also filter by level and melt
filter_data_rm_one = function(data, filter_parm, bad_val, level = NA, melt_flag = F)
{
  data2 = data
  data2 = data2[ data2[filter_parm] != bad_val & !is.na(data2[filter_parm]), ]
  
  if ( !is.na(level) )
  {
    level_taxa_pos = get_taxa_level_pos( data2, level )
    ids = names(data2);
    no_taxa_ids_pos = !grepl('__',ids) 
    # get only wanted level
    data2 = data2[,no_taxa_ids_pos | level_taxa_pos]
  }
  if (melt_flag)
  {
    ids = names(data2);
    data2 = melt(data2, id = ids[ !grepl('__',ids) ], variable_name = 'taxa', na.rm = F)
  }
  return(data2)
}


# get position of wanted level taxa in the data
get_taxa_level_pos = function( data, level )
{
  # take wanted level
  ids = names(data);
  
  # to get position of specific level get those who had x '__' but not x+1 '__' ('__' number in name matches the level)
  level_taxa_pos = grepl( sprintf('(__.*){%d}', level),ids ) & ! grepl( sprintf('(__.*){%d}', level+1),ids )
  
  return(level_taxa_pos)
}

# get position of wanted level taxa in the data. 
# filteres to metadta and wanted level (metadta optional by flag)
# relies on last metadta being "Description"
filter_taxa_to_level = function( data, level, metadata_flag = T )
{
  level_taxa_pos = get_taxa_level_pos( data, level )
  
  if (metadata_flag)
  {
    labs = names(taxa)
    pos = which(labs=='Description')
    return(data[,c(1:pos,which(level_taxa_pos))])
  }
  else
    return(data[,which(level_taxa_pos)])
}


# taking wanted level, ordering and melting data for plot and malting data 
set_taxa_for_plot = function( f_taxa, order_parm, taxa_level = 2, level_filter = T )
{
  # take wanted level
  ids = names(f_taxa);
  
  if ( level_filter ) {
    # to get position of specific level get those who had x '__' but not x+1 '__' ('__' number in name matches the level)
    level_taxa_pos = grepl( sprintf('(__.*){%d}', taxa_level),ids ) & ! grepl( sprintf('(__.*){%d}', taxa_level+1),ids )
  } else
    level_taxa_pos = grepl('__',ids) 
  # all bacteria name variables contain '__', take those who dont.
  no_taxa_ids_pos = !grepl('__',ids) 
  # get only wanted level
  plot_data = f_taxa[,no_taxa_ids_pos | level_taxa_pos]
  
  
  # melt dataframe into a form that ggplot can handel.
  
  plot_data = melt(plot_data, id = ids[no_taxa_ids_pos], variable_name = 'taxa', na.rm = F)
  # order by wanted parameter
  # plot_data = plot_data[order(plot_data[order_parm]), ]
  plot_data$SampleID = factor(plot_data$SampleID, levels = unique(plot_data$SampleID[ order(plot_data[[order_parm]]) ] ) )
  return(plot_data)
}


# ordering and melting data for plot
set_function_for_plot = function( f_taxa, order_parm, level_filter = F, taxa_level = 2 )
{
  # take wanted level
  ids = names(f_taxa);
  
  pos = which(ids == 'Description')
  level_taxa_pos = (pos+1):length(ids)
  # browser()
  if ( level_filter ) 
  {
    temp = which(check_id_level(ids) == taxa_level -1)
    level_taxa_pos = intersect(temp, level_taxa_pos) 
  }
  no_taxa_ids_pos = 1:pos
  # get only wanted level
  plot_data = f_taxa[,c(no_taxa_ids_pos, level_taxa_pos)]
  
  
  # melt dataframe into a form that ggplot can handel.
  plot_data = melt(plot_data, id = ids[no_taxa_ids_pos], variable_name = 'taxa', na.rm = F)
  
  # order by wanted parameter
  # plot_data = plot_data[order(plot_data[order_parm]), ]
  plot_data$SampleID = factor(plot_data$SampleID, levels = unique(plot_data$SampleID[ order(plot_data[[order_parm]]) ] ) )
  return(plot_data)
}

# get a string for title of area pragh contating all metadata filters. 
set_area_title = function(parameters, level = NA, x_prm, all_flag, average_flag = T)
{
  title = sprintf('Taxa L%d',level)
  
  title = sprintf('%s %s', title, parameters_str(parameters))
  title = sprintf('%sby %s', title, x_prm )
  if ( !all_flag )
  {
    if (average_flag)
      title = sprintf('%s Average', title)
    else
      title = sprintf('%s Median', title)
  }
  return(title)
}


# get a string for title contating all metadata filters. 
parameters_str = function(parameters)
{
  title = ''
  parm_names = names(parameters)
  for ( i in 1:length(parm_names) )
  {
    parm = unlist( parameters[parm_names[i]] )
    if ( length(parm) == 2 ) # in range, >= first and < second.
    {
      title = sprintf('%s%s=%d-%d ', title, parm_names[i], parm[1], parm[2] )
    }
    if ( length(parm) == 1 ) # equals
    {
      if ( !is.na(parm)  ) # can use NA for variables I dont want to check
      {
        if ( is.numeric(parm) )
          title = sprintf('%s %s=%d', title,parm_names[i], parm )
        else
          title = sprintf('%s %s=%s', title,parm_names[i], parm )
      }
    }
  }
  return(title)
}


# getting filtered data frame and ploting an area graph of it.
taxa_area_fig = function(taxa, x_parm ,tax_level = 2, x_precision = 1, all_flag = TRUE, filter_parameters = list())
{
  # ordering data for plot
  plot_data = set_taxa_for_plot( taxa, order_parm = x_parm, taxa_level = tax_level )
  
  # getting x labels
  lab = sort(taxa[[x_parm]])
  
  if (is.numeric( unlist( taxa[x_parm] ) ) )
  {
    lab = round(lab,x_precision )
  }
  # getting title
  ttl = set_area_title( filter_parameters, tax_level, x_parm, all_flag )
  # ploting
  if (all_flag)
  {
    p1 = ggplot(data = plot_data, aes(x = SampleID, y=value,group=taxa,fill=taxa)) +
      geom_area(position="fill", aes(colour = taxa, fill= taxa)) +
      scale_x_discrete(labels=lab , expand = c(0, 0)) +
      theme(axis.text = element_text(colour = 'black'), axis.title=element_text(size=18), plot.title = element_text(size=22)) +
      xlab(x_parm) + ylab('Taxa') + labs(title = ttl) +
      scale_y_continuous(expand = c(0, 0)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      theme(legend.position="bottom")
    # theme(legend.position="none")
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  }
  else
  {
    p1 = ggplot(data = plot_data, aes(y = value, x = plot_data[[ x_parm ]], colour = taxa, group = taxa)) + 
      stat_summary(fun.y = median, geom = "area", position="fill", aes(colour = taxa, fill= taxa)) +
      theme(axis.text = element_text(colour = 'black'), axis.title=element_text(size=18), plot.title = element_text(size=22)) +
      xlab(x_parm) + ylab('Taxa') + labs(title = ttl) + 
      theme(legend.position="bottom")
    if (is.numeric( unlist( f_taxa[x_parm] ) ) )
      p1 + scale_x_continuous(expand = c(0, 0)) 
    else
      p1 + scale_x_discrete(expand = c(0, 0)) 
  }
  
  return( list(p1,ttl, plot_data) )
}

# loading data, filtering and ploting an area graph by parametrs, all_flag is true is you want to plot all samples, for median false
plot_taxa_area = function(taxa_file, filter_parameters = list(), x_parm, tax_level = 2, na_strings = c('no_data','_','NA'), rm_na = TRUE, rm_parameters = NA, x_precision = 1, all_flag = TRUE, file_flag = TRUE, average_flag = T )
{
  ## loading data
  # for taxa read make sure to:
  # remove the # from header (from #sampleID)
  # replace all " " with "_"
  # file_flag = taxa_file is file to load if true, else it is a dataframe to use directly
  # average_flag = if True uses average(mean). if false uses median. only matters when all_flag is False
  if (file_flag)
    taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_strings)
  else
    taxa = taxa_file
  
  # remove taxa with low aboundences
  # browser()
  f_taxa = rmv_low_taxa(taxa, taxa_level = tax_level, cutoff = 0.01, high_num = 10)
  # filtering data
  if ( length(rm_met_parm)!=0 )
    if ( !is.na(rm_parameters) )
      f_taxa = filter_data_rm(f_taxa,rm_parameters, rmv_na = rm_na)
  f_taxa = filter_data(f_taxa,filter_parameters, rmv_na = rm_na)
  
  #removing NAs in tested parameter
  if (rm_na)
   f_taxa = filter_data_rm_one(f_taxa, x_parm, bad_val = '_')
  
  # ordering data for plot
  plot_data = set_taxa_for_plot( f_taxa, order_parm = x_parm, taxa_level = tax_level )
  # getting x labels
  lab = sort(f_taxa[[x_parm]])

  if (is.numeric( unlist( f_taxa[x_parm] ) ) )
  {
    lab = round(lab,x_precision )
  }
  # getting title
  ttl = set_area_title( filter_parameters, tax_level, x_parm, all_flag, average_flag )
  # ploting
  if (all_flag)
  {
    p1 = ggplot(data = plot_data, aes(x = SampleID, y=value, group=taxa, fill=taxa)) +
      geom_area(position="fill", aes(fill= taxa)) +
      scale_x_discrete(labels=lab , expand = c(0, 0)) +
      theme(axis.text = element_text(colour = 'black'), axis.title=element_text(size=18), plot.title = element_text(size=22)) +
      xlab(x_parm) + ylab('Taxa') + labs(title = ttl) +
      scale_y_continuous(expand = c(0, 0)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      theme(text = element_text(size=16), axis.text.x = element_text(size=14), 
            axis.text.y = element_text(size=14), axis.title = element_text(size = 16)) + 
      theme(legend.position="bottom") + #scale_fill_brewer(palette="Paired")
      scale_fill_manual(values=c("#660000", "#d53e4f", "#f46d43","#fdae61", 
                                 "#fee08b", "#ffffbf", "#e6f598","#abdda4", 
                                 "#66c2a5", "#3288bd", "#5e4fa2","#000033" )) 
    # theme(legend.position="none")
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  }
  else
  {
    p1 = ggplot(data = plot_data, aes(y = value, x = plot_data[[ x_parm ]], group = taxa)) 
    if (average_flag)
      p1 = p1 + stat_summary(fun.y = mean, geom = "area", position="fill", aes( fill= taxa))
    else
      p1 = p1 + stat_summary(fun.y = median, geom = "area", position="fill", aes( fill= taxa))
    p1 = p1 + theme(axis.text = element_text(colour = 'black'), axis.title=element_text(size=18), plot.title = element_text(size=22)) +
    xlab(x_parm) + ylab('Taxa') + labs(title = ttl) + 
    theme(text = element_text(size=16), axis.text.x = element_text(size=14), 
            axis.text.y = element_text(size=14), axis.title = element_text(size = 16)) + 
    theme(legend.position="bottom") + scale_y_continuous(expand = c(0, 0)) # +
    scale_fill_manual(values=c("#660000", "#d53e4f", "#f46d43","#fdae61", 
                               "#fee08b", "#ffffbf", "#e6f598","#abdda4", 
                               "#66c2a5", "#3288bd", "#5e4fa2","#000033" )) 
    if (is.numeric( unlist( f_taxa[x_parm] ) ) )
      p1 = p1 + scale_x_continuous(expand = c(0, 0)) 
    else
      p1 = p1 + scale_x_discrete(expand = c(0, 0)) 
  }
  
  return( list(p1,ttl, plot_data, f_taxa) )
}


# works for level 3 etc.
plot_taxa_area_lowLvl = function(taxa_file, filter_parameters = list(), x_parm, tax_level = 2, na_strings = c('no_data','_','NA'), rm_na = TRUE, rm_parameters = NA, x_precision = 1, all_flag = TRUE, file_flag = TRUE, average_flag = T )
{
  ## loading data
  # for taxa read make sure to:
  # remove the # from header (from #sampleID)
  # replace all " " with "_"
  # file_flag = taxa_file is file to load if true, else it is a dataframe to use directly
  # average_flag = if True uses average(mean). if false uses median. only matters when all_flag is False
  if (file_flag)
    taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_strings)
  else
    taxa = taxa_file
  
  # remove taxa with low aboundences
  # browser()
  f_taxa = rmv_low_taxa(taxa, taxa_level = tax_level, cutoff = 0.2, high_num = 20)
  # filtering data
  if ( length(rm_met_parm)!=0 )
    if ( !is.na(rm_parameters) )
      f_taxa = filter_data_rm(f_taxa,rm_parameters, rmv_na = rm_na)
  f_taxa = filter_data(f_taxa,filter_parameters, rmv_na = rm_na)
  
  #removing NAs in tested parameter
  if (rm_na)
    f_taxa = filter_data_rm_one(f_taxa, x_parm, bad_val = '_')
  
  # ordering data for plot
  plot_data = set_taxa_for_plot( f_taxa, order_parm = x_parm, taxa_level = tax_level )
  # getting x labels
  lab = sort(f_taxa[[x_parm]])
  
  if (is.numeric( unlist( f_taxa[x_parm] ) ) )
  {
    lab = round(lab,x_precision )
  }
  # getting title
  ttl = set_area_title( filter_parameters, tax_level, x_parm, all_flag, average_flag )
  # ploting
  if (all_flag)
  {
    p1 = ggplot(data = plot_data, aes(x = SampleID, y=value, group=taxa, fill=taxa)) +
      geom_area(position="fill", aes(fill= taxa)) +
      scale_x_discrete(labels=lab , expand = c(0, 0)) +
      theme(axis.text = element_text(colour = 'black'), axis.title=element_text(size=18), plot.title = element_text(size=22)) +
      xlab(x_parm) + ylab('Taxa') + labs(title = ttl) +
      scale_y_continuous(expand = c(0, 0)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      theme(text = element_text(size=16), axis.text.x = element_text(size=14), 
            axis.text.y = element_text(size=14), axis.title = element_text(size = 16)) + 
      theme(legend.position="bottom") # + #scale_fill_brewer(palette="Paired")
      # scale_fill_manual(values=c("#660000", "#d53e4f", "#f46d43","#fdae61", 
      #                            "#fee08b", "#ffffbf", "#e6f598","#abdda4", 
      #                            "#66c2a5", "#3288bd", "#5e4fa2","#000033" )) 
    # theme(legend.position="none")
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  }
  else
  {
    p1 = ggplot(data = plot_data, aes(y = value, x = plot_data[[ x_parm ]], group = taxa)) 
    if (average_flag)
      p1 = p1 + stat_summary(fun.y = mean, geom = "area", position="fill", aes( fill= taxa))
    else
      p1 = p1 + stat_summary(fun.y = median, geom = "area", position="fill", aes( fill= taxa))
    p1 = p1 + theme(axis.text = element_text(colour = 'black'), axis.title=element_text(size=18), plot.title = element_text(size=22)) +
      xlab(x_parm) + ylab('Taxa') + labs(title = ttl) + 
      theme(text = element_text(size=16), axis.text.x = element_text(size=14), 
            axis.text.y = element_text(size=14), axis.title = element_text(size = 16)) + 
      theme(legend.position="bottom") + scale_y_continuous(expand = c(0, 0)) # +
    # scale_fill_manual(values=c("#660000", "#d53e4f", "#f46d43","#fdae61", 
    #                            "#fee08b", "#ffffbf", "#e6f598","#abdda4", 
    #                            "#66c2a5", "#3288bd", "#5e4fa2","#000033" )) 
    if (is.numeric( unlist( f_taxa[x_parm] ) ) )
      p1 = p1 + scale_x_continuous(expand = c(0, 0)) 
    else
      p1 = p1 + scale_x_discrete(expand = c(0, 0)) 
  }
  
  return( list(p1,ttl, plot_data, f_taxa) )
}


set_taxa_for_area_plot = function(taxa_file, filter_parameters = list(), x_parm = NA, tax_level = 2, na_strings = c('no_data','_','NA'), rm_na = TRUE, rm_parameters = NA, x_precision = 1, all_flag = TRUE, file_flag = TRUE, average_flag = T, taxa_cutoff = 0.06, taxa_cutoff_high_num = 10, manual_color_flag = TRUE )
{
  ## loading data
  # for taxa read make sure to:
  # remove the # from header (from #sampleID)
  # replace all " " with "_"
  # file_flag = taxa_file is file to load if true, else it is a dataframe to use directly
  # average_flag = if True uses average(mean). if false uses median. only matters when all_flag is False
  if (file_flag)
    taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_strings)
  else
    taxa = taxa_file
  
  f_taxa = taxa
  # filtering data
  if ( length(rm_met_parm)!=0 )
    if ( !is.na(rm_parameters) )
      f_taxa = filter_data_rm(f_taxa,rm_parameters, rmv_na = rm_na)
  f_taxa = filter_data(f_taxa,filter_parameters, rmv_na = rm_na)
  
  #removing NAs in tested parameter
  if (rm_na)
    f_taxa = filter_data_rm_one(f_taxa, x_parm, bad_val = '_')
  
  # remove taxa with low aboundences
  # browser()
  f_taxa = rmv_low_taxa(f_taxa, taxa_level = tax_level, cutoff = taxa_cutoff, high_num = taxa_cutoff_high_num)
  # browser()
  
  # ordering data for plot
  plot_data = set_taxa_for_plot( f_taxa, order_parm = x_parm, taxa_level = tax_level )
  # getting x labels
  lab = sort(f_taxa[[x_parm]])
  
  if (is.numeric( unlist( f_taxa[x_parm] ) ) )
  {
    lab = round(lab,x_precision )
  }
  # getting title
  ttl = set_area_title( filter_parameters, tax_level, x_parm, all_flag, average_flag )
  # ploting
  if ( x_parm == 'Race' )
  {
    a<-levels(plot_data[[x_parm]])
    plot_data[[x_parm]] <- factor(plot_data[[x_parm]], levels = c(a[3],a[2],a[4],a[5], a[1]))
    
    names = c('Israeli_healthy' ,'Israeli_hospitalized','Malawi','Venezuela','America')
    new_names = c('Israeli healthy' ,'Israeli hospitalized','Malawi','Venezuela','American')
  }
  
  legs = levels(plot_data$taxa)
  new_legs = legs
  for (i in 1:length(legs))
    new_legs[i] = gsub('k__Bacteria.p__','',legs[i])
  
  return(list(plot_data, new_legs))
}

# loading data, filtering and ploting an area graph by parametrs, all_flag is true is you want to plot all samples, for median false
# v2 = new try with colors. specific for race median area plot
plot_taxa_area_v2 = function(taxa_file, filter_parameters = list(), x_parm, tax_level = 2, na_strings = c('no_data','_','NA'), rm_na = TRUE, rm_parameters = NA, x_precision = 1, all_flag = TRUE, file_flag = TRUE, average_flag = T, taxa_cutoff = 0.06, taxa_cutoff_high_num = 10, manual_color_flag = TRUE, order_parm = NA )
{
  ## loading data
  # for taxa read make sure to:
  # remove the # from header (from #sampleID)
  # replace all " " with "_"
  # file_flag = taxa_file is file to load if true, else it is a dataframe to use directly
  # average_flag = if True uses average(mean). if false uses median. only matters when all_flag is False
  if (file_flag)
    taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_strings)
  else
    taxa = taxa_file
  
  f_taxa = taxa
  # filtering data
  if ( length(rm_met_parm)!=0 )
    if ( !is.na(rm_parameters) )
      f_taxa = filter_data_rm(f_taxa,rm_parameters, rmv_na = rm_na)
  f_taxa = filter_data(f_taxa,filter_parameters, rmv_na = rm_na)
  
  #removing NAs in tested parameter
  if (rm_na)
    f_taxa = filter_data_rm_one(f_taxa, x_parm, bad_val = '_')
  
  # remove taxa with low aboundences
  # browser()
  f_taxa = rmv_low_taxa(f_taxa, taxa_level = tax_level, cutoff = taxa_cutoff, high_num = taxa_cutoff_high_num)
  # browser()
  
  # ordering data for plot
  plot_data = set_taxa_for_plot( f_taxa, order_parm = x_parm, taxa_level = tax_level )
  # getting x labels
  lab = sort(f_taxa[[x_parm]])
  
  if (is.numeric( unlist( f_taxa[x_parm] ) ) )
  {
    lab = round(lab,x_precision )
  }
  # getting title
  ttl = set_area_title( filter_parameters, tax_level, x_parm, all_flag, average_flag )
  # ploting
  if ( x_parm == 'Race' )
  {
    a<-levels(plot_data[[x_parm]])
    plot_data[[x_parm]] <- factor(plot_data[[x_parm]], levels = c(a[3],a[2],a[4],a[5], a[1]))
    
    names = c('Israeli_healthy' ,'Israeli_hospitalized','Malawi','Venezuela','America')
    new_names = c('Israeli healthy' ,'Israeli hospitalized','Malawi','Venezuela','American')
  }
  if ( x_parm == 'Cohort' )
  {
    a<-levels(plot_data[[x_parm]])
    plot_data[[x_parm]] <- factor(plot_data[[x_parm]], levels = c(a[5],a[4],a[3],a[6], a[7],a[1],a[2]))
    
    names = c('Israeli_Hospitalized','Israeli_Healthy','Israeli_Control','Malawi_Healthy','Venezuela_Healthy','America_Healthy','America_healthy_students')
    new_names = c('Israeli hospitalized','Israeli healthy 1','Israeli healthy 2','Malawians','Amerindians','US' ,'US students')
  }
  
  legs = levels(plot_data$taxa)
  new_legs = legs
  for (i in 1:length(legs))
    new_legs[i] = gsub('k__Bacteria.p__','',legs[i])
  
  if (all_flag)
  {
    p1 = ggplot(data = plot_data, aes(x = SampleID, y=value, group=taxa, fill=taxa)) +
      geom_area(position="fill", aes(fill= taxa)) +
      scale_x_discrete(labels=lab , expand = c(0, 0)) +
      theme(axis.text = element_text(colour = 'black'), axis.title=element_text(size=18), plot.title = element_text(size=22)) +
      xlab(x_parm) + ylab('Taxa') + labs(title = ttl) +
      scale_y_continuous(expand = c(0, 0)) + 
      scale_fill_discrete( labels=new_legs, breaks = legs, name = 'Phylum')
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      theme(text = element_text(size=16), axis.text.x = element_text(size=14), 
            axis.text.y = element_text(size=14), axis.title = element_text(size = 16)) + 
      theme(legend.position="bottom") # + #  scale_fill_brewer(palette="Paired")
      # scale_fill_manual(values=c("#660000", "#d53e4f", "#f46d43","#fdae61", 
      #                           "#fee08b", "#ffffbf", "#e6f598","#abdda4", 
      #                           "#66c2a5", "#3288bd", "#5e4fa2","#000033" )) 
    # theme(legend.position="none")
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  }
  else
  {
    p1 = ggplot(data = plot_data, aes(y = value, x = plot_data[[ x_parm ]], group = taxa)) 
    if (average_flag)
      p1 = p1 + stat_summary(fun.y = mean, geom = "area", position="fill", aes( fill= taxa))
    else
      # p1 = p1 + stat_summary(fun.y = median, geom = "bar", position="fill", aes( fill= taxa))
      p1 = p1 + stat_summary(fun.y = median, geom = "area", position="fill", aes( fill= taxa))
    p1 = p1 + theme(axis.text = element_text(colour = 'black'), axis.title=element_text(size=18), plot.title = element_text(size=22)) +
      xlab(x_parm) + ylab('Taxa') + labs(title = ttl) + 
      theme(text = element_text(size=14), axis.text.x = element_text(size=13), 
            axis.text.y = element_text(size=13), axis.title = element_text(size = 16)) + 
      theme(legend.position="bottom") + scale_y_continuous(expand = c(0, 0)) + theme_classic() 
      # scale_fill_manual(values=  color_gradiant( c("blue","red"),1:length(levels(plot_data$taxa)) ) ) + 
      if ( manual_color_flag ) {
        p1 = p1 + scale_fill_manual( values = c("#000033","#3288bd","#66c2a5", "#fee08b", "#f46d43", "#d53e4f","#660000"),
                         labels=new_legs, breaks = legs, name = 'Phylum' ) 
      } else {
        p1 = p1 + scale_fill_discrete( labels=new_legs, breaks = legs, name = 'Phylum')
      }
    
      #scale_fill_manual(values=c("#660000", "#d53e4f", "#f46d43","#fdae61", 
      #                           "#fee08b", "#ffffbf", "#e6f598","#abdda4")) #, 
      #                          "#66c2a5", "#3288bd", "#5e4fa2","#000033" )) 
    if (is.numeric( unlist( f_taxa[x_parm] ) ) )
      p1 = p1 + scale_x_continuous(expand = c(0, 0)) 
    else
    {
      if ( x_parm == 'Race' | x_parm == 'Cohort' )
        p1 = p1 + scale_x_discrete(expand = c(0, 0), breaks=names, labels=new_names)
      else
        p1 = p1 + scale_x_discrete(expand = c(0, 0)) 
    }
  }
  
  return( list(p1,ttl, plot_data, f_taxa) )
}

# loading data, filtering and ploting an area graph by parametrs, all_flag is true is you want to plot all samples, for median false
# v3 = added labels dataframe optianal input for lables.
plot_taxa_area_v3 = function(taxa_file, filter_parameters = list(), x_parm, tax_level = 2, na_strings = c('no_data','_','NA','na'), rm_na = TRUE, rm_parameters = NA, x_precision = 1, all_flag = TRUE, file_flag = TRUE, average_flag = F, taxa_cutoff = 0.06, taxa_cutoff_high_num = 10, manual_color_flag = TRUE, order_parm = NA, labelss = NA, cols = c("#000033","#3288bd","#66c2a5", "#fee08b", "#f46d43", "#d53e4f","#660000") )
{
  ## loading data
  # for taxa read make sure to:
  # remove the # from header (from #sampleID)
  # replace all " " with "_"
  # file_flag = taxa_file is file to load if true, else it is a dataframe to use directly
  # average_flag = if True uses average(mean). if false uses median. only matters when all_flag is False
  if (file_flag)
    taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_strings)
  else
    taxa = taxa_file
  
  f_taxa = taxa
  # filtering data
  if ( length(rm_met_parm)!=0 )
    if ( !is.na(rm_parameters) )
      f_taxa = filter_data_rm(f_taxa,rm_parameters, rmv_na = rm_na)
  f_taxa = filter_data(f_taxa,filter_parameters, rmv_na = rm_na)
  
  #removing NAs in tested parameter
  if (rm_na)
    f_taxa = filter_data_rm_one(f_taxa, x_parm, bad_val = '_')
  
  # remove taxa with low aboundences
  # browser()
  f_taxa = rmv_low_taxa(f_taxa, taxa_level = tax_level, cutoff = taxa_cutoff, high_num = taxa_cutoff_high_num)
  # browser()
  
  # ordering data for plot
  plot_data = set_taxa_for_plot( f_taxa, order_parm = x_parm, taxa_level = tax_level )
  # getting x labels
  lab = sort(f_taxa[[x_parm]])
  
  
  if ( is.numeric( unlist( f_taxa[x_parm] ) ) )
  {
    lab = round(lab,x_precision )
  }
  # getting title
  ttl = set_area_title( filter_parameters, tax_level, x_parm, all_flag, average_flag )
  # ploting
  
  # if want to change names use a labels dataframe with all the neccecary data
  if ( !is.na(labelss) )
  {
    if ( check_val == labelss$check_val )
    {
      a<-levels(plot_data[[x_parm]])
      # browser()
      lvls = c()
      for (ord in labelss$names_order )
        lvls =  c(lvls, a[ord])
      plot_data[[x_parm]] <- factor(plot_data[[x_parm]], levels = lvls)
      
      names = labelss$names
      new_names = labelss$new_names
    }
  }
  
  legs = levels(plot_data$taxa)
  new_legs = legs
  for (i in 1:length(legs))
    new_legs[i] = gsub('k__Bacteria.p__','',legs[i])
  
  if (all_flag)
  {
    p1 = ggplot(data = plot_data, aes(x = SampleID, y=value, group=taxa, fill=taxa)) +
      geom_area(position="fill", aes(fill= taxa)) +
      scale_x_discrete(labels=lab , expand = c(0, 0)) +
      theme(axis.text = element_text(colour = 'black'), axis.title=element_text(size=18), plot.title = element_text(size=22)) +
      xlab(x_parm) + ylab('Taxa') + labs(title = ttl) +
      scale_y_continuous(expand = c(0, 0)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      theme(text = element_text(size=16), axis.text.x = element_text(size=14), 
            axis.text.y = element_text(size=14), axis.title = element_text(size = 16)) + 
      scale_fill_discrete( labels=new_legs, breaks = legs, name = 'Phylum' ) +
      theme(legend.position="bottom") # + #  scale_fill_brewer(palette="Paired")
    # scale_fill_manual(values=c("#660000", "#d53e4f", "#f46d43","#fdae61", 
    #                           "#fee08b", "#ffffbf", "#e6f598","#abdda4", 
    #                           "#66c2a5", "#3288bd", "#5e4fa2","#000033" )) 
    # theme(legend.position="none")
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  }
  else
  {
    p1 = ggplot(data = plot_data, aes(y = value, x = plot_data[[ x_parm ]], group = taxa)) 
    if (average_flag)
      p1 = p1 + stat_summary(fun.y = mean, geom = "area", position="fill", aes( fill= taxa))
    else
      # p1 = p1 + stat_summary(fun.y = median, geom = "bar", position="fill", aes( fill= taxa))
      p1 = p1 + stat_summary(fun.y = median, geom = "area", position="fill", aes( fill= taxa))
    p1 = p1 + theme(axis.text = element_text(colour = 'black'), axis.title=element_text(size=18), plot.title = element_text(size=22)) +
      xlab(x_parm) + ylab('Taxa') + labs(title = ttl) + 
      theme(text = element_text(size=14), axis.text.x = element_text(size=13), 
            axis.text.y = element_text(size=13), axis.title = element_text(size = 16)) + 
      theme(legend.position="bottom") + scale_y_continuous(expand = c(0, 0)) + theme_classic() 
    # scale_fill_manual(values=  color_gradiant( c("blue","red"),1:length(levels(plot_data$taxa)) ) ) + 
    if ( manual_color_flag ) {
      p1 = p1 + scale_fill_manual( values = cols,
                                   labels=new_legs, breaks = legs, name = 'Phylum' ) 
    } else {
      p1 = p1 + scale_fill_discrete( labels=new_legs, breaks = legs, name = 'Phylum')
    }
    
    #scale_fill_manual(values=c("#660000", "#d53e4f", "#f46d43","#fdae61", 
    #                           "#fee08b", "#ffffbf", "#e6f598","#abdda4")) #, 
    #                          "#66c2a5", "#3288bd", "#5e4fa2","#000033" )) 
    if (is.numeric( unlist( f_taxa[x_parm] ) ) )
      p1 = p1 + scale_x_continuous(expand = c(0, 0)) 
    else
    {
      p1 = p1 + scale_x_discrete(expand = c(0, 0)) 
      if ( !is.na(labelss)  )
        if ( check_val == labelss$check_val )
          p1 = p1 + scale_x_discrete(expand = c(0, 0), breaks=names, labels=new_names)

        
    }
  }
  
  return( list(p1,ttl, plot_data, f_taxa) )
}

# loading data, filtering and ploting an area graph by parametrs, all_flag is true is you want to plot all samples, for median false
# v3 = added labels dataframe optianal input for lables.
plot_taxa_bar = function(taxa_file, filter_parameters = list(), x_parm, tax_level = 2, na_strings = c('no_data','_','NA','na'), rm_na = TRUE, rm_parameters = NA, x_precision = 1, all_flag = TRUE, file_flag = TRUE, average_flag = F, taxa_cutoff = 0.06, taxa_cutoff_high_num = 10, manual_color_flag = TRUE, order_parm = NA, labelss = NA, cols = c("#000033","#3288bd","#66c2a5", "#fee08b", "#f46d43", "#d53e4f","#660000") )
{
  ## loading data
  # for taxa read make sure to:
  # remove the # from header (from #sampleID)
  # replace all " " with "_"
  # file_flag = taxa_file is file to load if true, else it is a dataframe to use directly
  # average_flag = if True uses average(mean). if false uses median. only matters when all_flag is False
  if (file_flag)
    taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_strings)
  else
    taxa = taxa_file
  f_taxa = taxa
  # filtering data
  if ( length(rm_met_parm)!=0 )
    if ( !is.na(rm_parameters) )
      f_taxa = filter_data_rm(f_taxa,rm_parameters, rmv_na = rm_na)
  f_taxa = filter_data(f_taxa,filter_parameters, rmv_na = rm_na)
  
  #removing NAs in tested parameter
  if (rm_na)
    f_taxa = filter_data_rm_one(f_taxa, x_parm, bad_val = '_')
  
  # remove taxa with low aboundences
  # browser()
  f_taxa = rmv_low_taxa(f_taxa, taxa_level = tax_level, cutoff = taxa_cutoff, high_num = taxa_cutoff_high_num)
  # browser()
  
  # ordering data for plot
  plot_data = set_taxa_for_plot( f_taxa, order_parm = x_parm, taxa_level = tax_level )
  # getting x labels
  lab = sort(f_taxa[[x_parm]])
  
  
  if ( is.numeric( unlist( f_taxa[x_parm] ) ) )
  {
    lab = round(lab,x_precision )
  }
  # getting title
  ttl = set_area_title( filter_parameters, tax_level, x_parm, all_flag, average_flag )
  # ploting
  
  # if want to change names use a labels dataframe with all the neccecary data
  if ( !is.na(labelss) )
  {
    if ( check_val == labelss$check_val )
    {
      a<-levels(plot_data[[x_parm]])
      # browser()
      lvls = c()
      for (ord in labelss$names_order )
        lvls =  c(lvls, a[ord])
      plot_data[[x_parm]] <- factor(plot_data[[x_parm]], levels = lvls)
      
      names = labelss$names
      new_names = labelss$new_names
    }
  }
  
  legs = levels(plot_data$taxa)
  new_legs = legs
  for (i in 1:length(legs))
    new_legs[i] = gsub('k__Bacteria.p__','',legs[i])
  
  if (all_flag)
  {
    p1 = ggplot(data = plot_data, aes(x = SampleID, y=value, group=taxa, fill=taxa)) +
      geom_bar(position="fill", aes(fill= taxa)) +
      scale_x_discrete(labels=lab , expand = c(0, 0)) +
      theme(axis.text = element_text(colour = 'black'), axis.title=element_text(size=18), plot.title = element_text(size=22)) +
      xlab(x_parm) + ylab('Taxa') + labs(title = ttl) +
      scale_y_continuous(expand = c(0, 0)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      theme(text = element_text(size=16), axis.text.x = element_text(size=14), 
            axis.text.y = element_text(size=14), axis.title = element_text(size = 16)) + 
      scale_fill_discrete( labels=new_legs, breaks = legs, name = 'Phylum' ) +
      theme(legend.position="bottom") # + #  scale_fill_brewer(palette="Paired")
    # scale_fill_manual(values=c("#660000", "#d53e4f", "#f46d43","#fdae61", 
    #                           "#fee08b", "#ffffbf", "#e6f598","#abdda4", 
    #                           "#66c2a5", "#3288bd", "#5e4fa2","#000033" )) 
    # theme(legend.position="none")
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  }
  else
  {
    p1 = ggplot(data = plot_data, aes(y = value, x = plot_data[[ x_parm ]], group = taxa)) 
    if (average_flag)
      p1 = p1 + stat_summary(fun.y = mean, geom = "area", position="fill", aes( fill= taxa))
    else
      # p1 = p1 + stat_summary(fun.y = median, geom = "bar", position="fill", aes( fill= taxa))
      p1 = p1 + stat_summary(fun.y = median, geom = "bar", position="fill", aes( fill= taxa))
    p1 = p1 + theme(axis.text = element_text(colour = 'black'), axis.title=element_text(size=18), plot.title = element_text(size=22)) +
      xlab(x_parm) + ylab('Taxa') + labs(title = ttl) + 
      theme(text = element_text(size=14), axis.text.x = element_text(size=13), 
            axis.text.y = element_text(size=13), axis.title = element_text(size = 16)) + 
      theme(legend.position="bottom") + scale_y_continuous(expand = c(0, 0)) + theme_classic() 
    # scale_fill_manual(values=  color_gradiant( c("blue","red"),1:length(levels(plot_data$taxa)) ) ) + 
    if ( manual_color_flag ) {
      p1 = p1 + scale_fill_manual( values = cols,
                                   labels=new_legs, breaks = legs, name = 'Phylum' ) 
    } else {
      p1 = p1 + scale_fill_discrete( labels=new_legs, breaks = legs, name = 'Phylum')
    }
    
    #scale_fill_manual(values=c("#660000", "#d53e4f", "#f46d43","#fdae61", 
    #                           "#fee08b", "#ffffbf", "#e6f598","#abdda4")) #, 
    #                          "#66c2a5", "#3288bd", "#5e4fa2","#000033" )) 
    if (is.numeric( unlist( f_taxa[x_parm] ) ) )
      p1 = p1 + scale_x_continuous(expand = c(0, 0)) 
    else
    {
      p1 = p1 + scale_x_discrete(expand = c(0, 0)) 
      if ( !is.na(labelss)  )
        if ( check_val == labelss$check_val )
          p1 = p1 + scale_x_discrete(expand = c(0, 0), breaks=names, labels=new_names)
        
        
    }
  }
  
  return( list(p1,ttl, plot_data, f_taxa) )
}

# v2 = rewrite.
plot_taxa_bar_v2 = function(taxa_file, x_parm, tax_level = 2, na_strings = c('no_data','_','NA','na'), rm_na = TRUE, rm_parameters = NA, x_precision = 1, all_flag = TRUE, file_flag = TRUE, average_flag = F, taxa_cutoff = 0.06, taxa_cutoff_high_num = 10, manual_color_flag = TRUE, order_parm = NA, labelss = NA)
{
  ## loading data
  # for taxa read make sure to:
  # remove the # from header (from #sampleID)
  # replace all " " with "_"
  # file_flag = taxa_file is file to load if true, else it is a dataframe to use directly
  # average_flag = if True uses average(mean). if false uses median. only matters when all_flag is False
  if (file_flag)
    taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_strings)
  else
    taxa = taxa_file
  
  f_taxa = taxa
  # filtering data
  
  #removing NAs in tested parameter
  if (rm_na)
    f_taxa = filter_data_rm_one(f_taxa, x_parm, bad_val = '_')
  
  # remove taxa with low aboundences
  f_taxa = rmv_low_taxa( f_taxa, taxa_level = tax_level, cutoff = taxa_cutoff, high_num = taxa_cutoff_high_num )
  
  # ordering data for plot
  plot_data = set_taxa_for_plot( f_taxa, order_parm = x_parm, taxa_level = tax_level )
  # getting x labels
  lab = sort(f_taxa[[x_parm]])
  
  if ( is.numeric( unlist( f_taxa[x_parm] ) ) )
  {
    lab = round(lab,x_precision )
  }
  # getting title
  ttl = set_area_title(  list() , tax_level, x_parm, all_flag, average_flag )
  # ploting
  
  # if want to change names use a labels dataframe with all the neccecary data
  if ( !is.na(labelss) )
  {
    if ( check_val == labelss$check_val )
    {
      a<-levels(plot_data[[x_parm]])
      # browser()
      lvls = c()
      for (ord in labelss$names_order )
        lvls =  c(lvls, a[ord])
      plot_data[[x_parm]] <- factor(plot_data[[x_parm]], levels = lvls)
      
      names = labelss$names
      new_names = labelss$new_names
    }
  }
  
  legs = levels(plot_data$taxa)
  new_legs = legs
  for (i in 1:length(legs))
    new_legs[i] = gsub('k__Bacteria.p__','',legs[i])
  if (all_flag)
  {
    p1 = ggplot(data = plot_data, aes(x = SampleID, y=value, group=taxa, fill=taxa)) +
      geom_bar(position="fill", aes(fill= taxa)) +
      scale_x_discrete(labels=lab , expand = c(0, 0)) +
      theme(axis.text = element_text(colour = 'black'), axis.title=element_text(size=18), plot.title = element_text(size=22)) +
      xlab(x_parm) + ylab('Taxa') + labs(title = ttl) +
      scale_y_continuous(expand = c(0, 0)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      theme(text = element_text(size=16), axis.text.x = element_text(size=14), 
            axis.text.y = element_text(size=14), axis.title = element_text(size = 16)) + 
      scale_fill_discrete( labels=new_legs, breaks = legs, name = 'Phylum' ) +
      theme(legend.position="bottom") # + #  scale_fill_brewer(palette="Paired")
    # scale_fill_manual(values=c("#660000", "#d53e4f", "#f46d43","#fdae61", 
    #                           "#fee08b", "#ffffbf", "#e6f598","#abdda4", 
    #                           "#66c2a5", "#3288bd", "#5e4fa2","#000033" )) 
    # theme(legend.position="none")
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  }
  else
  {
    p1 = ggplot(data = plot_data, aes(y = value, x = plot_data[[ x_parm ]], group = taxa)) 
    if (average_flag)
      p1 = p1 + stat_summary(fun.y = mean, geom = "area", position="fill", aes( fill= taxa))
    else
      # p1 = p1 + stat_summary(fun.y = median, geom = "bar", position="fill", aes( fill= taxa))
      p1 = p1 + stat_summary(fun.y = median, geom = "bar", position="fill", aes( fill= taxa))
    p1 = p1 + theme(axis.text = element_text(colour = 'black'), axis.title=element_text(size=18), plot.title = element_text(size=22)) +
      xlab(x_parm) + ylab('Taxa') + labs(title = ttl) + 
      theme(text = element_text(size=14), axis.text.x = element_text(size=13), 
            axis.text.y = element_text(size=13), axis.title = element_text(size = 16)) + 
      theme(legend.position="bottom") + scale_y_continuous(expand = c(0, 0)) + theme_classic() 
    # scale_fill_manual(values=  color_gradiant( c("blue","red"),1:length(levels(plot_data$taxa)) ) ) + 
    if ( manual_color_flag ) {
      # p1 = p1 + scale_fill_manual( values = cols,
      #                              labels=new_legs, breaks = legs, name = 'Phylum' ) 
      p1 = p1 + scale_fill_brewer( palette = "Spectral",
                                    labels=new_legs, breaks = legs, name = 'Phylum' ) 
    } else {
      p1 = p1 + scale_fill_discrete( labels=new_legs, breaks = legs, name = 'Phylum')
    }
    
    #scale_fill_manual(values=c("#660000", "#d53e4f", "#f46d43","#fdae61", 
    #                           "#fee08b", "#ffffbf", "#e6f598","#abdda4")) #, 
    #                          "#66c2a5", "#3288bd", "#5e4fa2","#000033" )) 
    if (is.numeric( unlist( f_taxa[x_parm] ) ) )
      p1 = p1 + scale_x_continuous(expand = c(0, 0)) 
    else
    {
      p1 = p1 + scale_x_discrete(expand = c(0, 0)) 
      if ( !is.na(labelss)  )
        if ( check_val == labelss$check_val )
          p1 = p1 + scale_x_discrete(expand = c(0, 0), breaks=names, labels=new_names)
    }
  }
  
  return( list(p1,ttl, plot_data, f_taxa) )
}
# loading data, filtering and ploting an area graph by parametrs, all_flag is true is you want to plot all samples, for median false
# v3 = added labels dataframe optianal input for lables.
plot_taxa_pie = function(taxa, tax_level = 2, na_strings = c('no_data','_','NA','na'), taxa_cutoff = 0.06, taxa_cutoff_high_num = 10, all_flag =F )
{
  x_precision = 1

  f_taxa = taxa
  
  x_parm = 'SampleID'
  # ordering data for plot
  plot_data = set_taxa_for_plot( f_taxa, order_parm = x_parm, taxa_level = tax_level )
  # getting x labels
  lab = sort(f_taxa[[x_parm]])
  
  # getting title
  # ttl = set_area_title( filter_parameters, tax_level, x_parm, all_flag, average_flag )
  # ploting
  
  plot_data$value=plot_data$value*100
  
  legs = levels(plot_data$taxa)
  new_legs = legs
  for (i in 1:length(legs))
    new_legs[i] = gsub('k__Bacteria.p__','',legs[i])
  
  p1 = ggplot(plot_data, aes(x="", y=value, fill=taxa))+
   geom_bar(width = 1, stat = "identity") +  coord_polar("y", start = 0)
  

  p1 = p1 + theme(axis.text = element_text(colour = 'black'), axis.title=element_text(size=18), plot.title = element_text(size=22)) +
    xlab('') + ylab('') + # labs(title = ttl) + theme(legend.position="bottom") +
      theme_minimal() 
  p1 = p1 + scale_fill_brewer( labels=new_legs, breaks = legs, name = 'Phylum',palette="Spectral")
  
  
  return( list(p1, plot_data, f_taxa) )
}

plot_taxa_bar_v3 = function(taxa, x_parm, tax_level = 2, na_strings = c('no_data','_','NA','na'), average_flag = F, manual_color_flag = TRUE, cols = c("#000033","#3288bd","#66c2a5", "#fee08b", "#f46d43", "#d53e4f","#660000"), high_taxa_num = 6,all_flag = F, x_precision = 2, as_character_flag =T )
{
  f_taxa = taxa
  if (as_character_flag)
  {
    f_taxa[[x_parm]] = as.character(f_taxa[[x_parm]])
    # adding n= to x labs
    temp = f_taxa[[x_parm]]
    for ( i in 1:length(f_taxa[[x_parm]] ) )
    {
      f_taxa[[x_parm]][i] = sprintf('%s (n=%d)', temp[i], sum( temp == temp[i], na.rm = T ) )
      if (is.na(temp[i]))
        f_taxa[[x_parm]][i] = sprintf('%s (n=%d)', temp[i], sum( is.na(temp) ) )
    }
  }
  # changing taxa with low abundances to "other" (after high_taxa_num most abundant bacteria )
  f_taxa = rmv_low_taxa_v2(f_taxa, taxa_level = tax_level, high_taxa_num = high_taxa_num)
  # ordering data for plot
  plot_data = set_taxa_for_plot( f_taxa, order_parm = x_parm, taxa_level = tax_level, level_filter = F )
  # getting x labels
  lab = sort(f_taxa[[x_parm]])
  if ( is.numeric( unlist( f_taxa[x_parm] ) ) )
    lab = round(lab,x_precision )
  ttl = set_area_title( parameters = list(), tax_level, x_parm, all_flag, average_flag )
  legs = levels(plot_data$taxa)
  new_legs = legs
  for (i in 1:length(legs))
  {
    new_legs[i] = gsub('k__Bacteria.p__','',legs[i])
    new_legs[i] = gsub('k__','',new_legs[i])
  }
  if (all_flag)
  {
    p1 = ggplot(data = plot_data, aes(x = SampleID, y=value, group=taxa, fill=taxa)) +
      geom_bar(position="fill", aes(fill= taxa), stat='identity') +
      scale_x_discrete(labels=lab , expand = c(0, 0)) +
      theme(axis.text = element_text(colour = 'black'),
            axis.title=element_text(size=18), plot.title = element_text(size=22)) +
      xlab(x_parm) + ylab('Taxa') + labs(title = ttl) +
      scale_y_continuous(expand = c(0, 0)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      theme(text = element_text(size=16), axis.text.x = element_text(size=14), 
            axis.text.y = element_text(size=14), axis.title = element_text(size = 16)) + 
      scale_fill_manual( values = cols, labels=new_legs,
                         breaks = legs, name = 'Phylum' ) 
      # scale_fill_discrete( labels=new_legs, breaks = legs, name = 'Phylum' ) # +
      # theme(legend.position="bottom") # + #  scale_fill_brewer(palette="Paired")
    # scale_fill_manual(values=c("#660000", "#d53e4f", "#f46d43","#fdae61", 
    #                           "#fee08b", "#ffffbf", "#e6f598","#abdda4", 
    #                           "#66c2a5", "#3288bd", "#5e4fa2","#000033" )) 
    # theme(legend.position="none")
  }
  else
  {
    p1 = ggplot(data = plot_data, aes_string(y = 'value', x =x_parm, group = 'taxa')) 
    if (average_flag)
    {
      p1 = p1 + stat_summary(fun.y = mean, geom = "bar", position="fill", aes( fill= taxa))
    } else
      # p1 = p1 + stat_summary(fun.y = median, geom = "bar", position="fill", aes( fill= taxa))
      p1 = p1 + stat_summary(fun = median, geom = "bar", position="fill", aes( fill= taxa))
    p1 = p1 + theme(axis.text = element_text(colour = 'black'), axis.title=element_text(size=18), plot.title = element_text(size=22)) +
      xlab(x_parm) + ylab('Taxa') + labs(title = ttl) + 
      scale_y_continuous(expand = c(0, 0)) + 
      theme_classic() + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    # scale_fill_manual(values=  color_gradiant( c("blue","red"),1:length(levels(plot_data$taxa)) ) ) + 
    if ( manual_color_flag ) {
      p1 = p1 + scale_fill_manual( values = cols,
                                   labels=new_legs, breaks = legs, name = 'Phylum' ) 
    } else {
      p1 = p1 + scale_fill_discrete( labels=new_legs, breaks = legs, name = 'Phylum')
    }
    
    #scale_fill_manual(values=c("#660000", "#d53e4f", "#f46d43","#fdae61", 
    #                           "#fee08b", "#ffffbf", "#e6f598","#abdda4")) #, 
    #                          "#66c2a5", "#3288bd", "#5e4fa2","#000033" )) 
    if (is.numeric( unlist( f_taxa[x_parm] ) ) )
      p1 = p1 + scale_x_continuous(expand = c(0, 0)) 
    else
    {
      p1 = p1 + scale_x_discrete(expand = c(0, 0)) 
    }
  }
  return( list(p1,ttl, plot_data, f_taxa) )
}



# loading data, filtering and ploting an area graph by parametrs, all_flag is true is you want to plot all samples, for median false
# v2 = new try with colors. specific for race median area plot
plot_picrust_area = function(taxa_file, filter_parameters, x_parm, na_strings = c('no_data','_','NA'), rm_na = TRUE, rm_parameters = NA, x_precision = 1, all_flag = TRUE, file_flag = TRUE, average_flag = T, taxa_cutoff = 0.06, taxa_cutoff_high_num = 10, manual_color_flag = TRUE, tax_level=2 )
{
  ## loading data
  # for taxa read make sure to:
  # remove the # from header (from #sampleID)
  # replace all " " with "_"
  # file_flag = taxa_file is file to load if true, else it is a dataframe to use directly
  # average_flag = if True uses average(mean). if false uses median. only matters when all_flag is False
  if (file_flag)
    taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_strings)
  else
    taxa = taxa_file
  
  f_taxa = taxa
  # filtering data
  if ( length(rm_met_parm)!=0 )
    if ( !is.na(rm_parameters) )
      f_taxa = filter_data_rm(f_taxa,rm_parameters, rmv_na = rm_na)
  f_taxa = filter_data(f_taxa,filter_parameters, rmv_na = rm_na)
  
  #removing NAs in tested parameter
  if (rm_na)
    f_taxa = filter_data_rm_one(f_taxa, x_parm, bad_val = '_')
  
  # remove taxa with low aboundences
  # browser()
  # f_taxa = rmv_low_function(f_taxa, cutoff = taxa_cutoff, high_num = taxa_cutoff_high_num)
  # browser()
  f_taxa = f_taxa[!is.na(f_taxa$SampleID),]
  
  # ordering data for plot
  plot_data = set_function_for_plot( f_taxa, order_parm = x_parm, level_filter = T, taxa_level = tax_level )
  # getting x labels
  lab = sort(f_taxa[[x_parm]])
  
  if (is.numeric( unlist( f_taxa[x_parm] ) ) )
  {
    lab = round(lab,x_precision )
  }
  # getting title
  ttl = set_area_title( filter_parameters, tax_level, x_parm, all_flag, average_flag )
  ttl = gsub('Taxa','Metagenomic Function',ttl)
  # ttl = sprintf('Metagenomic Function L%d by %s', tax_level, x_parm)
  # if ( !all_flag )
  #   if ( average_flag ){
  #     ttl = sprintf('%s Average', ttl)
  #   } else {
  #     ttl = sprintf('%s Median', ttl)
  #   }
  
  # ploting
  if ( x_parm == 'Race' )
  {
    a<-levels(plot_data[[x_parm]])
    plot_data[[x_parm]] <- factor(plot_data[[x_parm]], levels = c(a[3],a[2],a[4],a[5], a[1]))
    
    names = c('Israeli_healthy' ,'Israeli_hospitalized','Malawi','Venezuela','America')
    new_names = c('Israeli healthy' ,'Israeli hospitalized','Malawi','Venezuela','American')
  }
  
  legs = levels(plot_data$taxa)
  new_legs = legs
  for (i in 1:length(legs))
    new_legs[i] = gsub('k__Bacteria.p__','',legs[i])
  if (all_flag)
  {
    p1 = ggplot(data = plot_data, aes(x = SampleID, y=value, group=taxa, fill=taxa)) +
      geom_area(position="fill", aes(fill= taxa)) +
      scale_x_discrete(labels=lab , expand = c(0, 0)) +
      theme(axis.text = element_text(colour = 'black'), axis.title=element_text(size=18), plot.title = element_text(size=22)) +
      xlab(x_parm) + ylab('Abundance') + labs(title = ttl) +
      scale_y_continuous(expand = c(0, 0)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      theme(text = element_text(size=16), axis.text.x = element_text(size=14), 
            axis.text.y = element_text(size=14), axis.title = element_text(size = 16)) + 
      theme(legend.position="right") # + #  scale_fill_brewer(palette="Paired")
    if ( manual_color_flag ) {
      p1 = p1 + scale_fill_manual( values = c("#000033","#3288bd","#66c2a5", "#fee08b", "#f46d43", "#d53e4f","#660000"),
                                   labels=new_legs, breaks = legs, name = 'Function' ) 
    } else {
      p1 = p1 + scale_fill_discrete( labels=new_legs, breaks = legs, name = 'Function')
    }
    # scale_fill_manual(values=c("#660000", "#d53e4f", "#f46d43","#fdae61", 
    #                           "#fee08b", "#ffffbf", "#e6f598","#abdda4", 
    #                           "#66c2a5", "#3288bd", "#5e4fa2","#000033" )) 
    # theme(legend.position="none")
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  }
  else
  {
    p1 = ggplot(data = plot_data, aes(y = value, x = plot_data[[ x_parm ]], group = taxa)) 
    if (average_flag)
      p1 = p1 + stat_summary(fun.y = mean, geom = "area", position="fill", aes( fill= taxa))
    else
      # p1 = p1 + stat_summary(fun.y = median, geom = "bar", position="fill", aes( fill= taxa))
      p1 = p1 + stat_summary(fun.y = median, geom = "area", position="fill", aes( fill= taxa))
    p1 = p1 + theme(axis.text = element_text(colour = 'black'), axis.title=element_text(size=18), plot.title = element_text(size=22)) +
      xlab(x_parm) + ylab('Abundance') + labs(title = ttl) + 
      theme(text = element_text(size=14), axis.text.x = element_text(size=13), 
            axis.text.y = element_text(size=13), axis.title = element_text(size = 16)) + 
      theme(legend.position="bottom") + scale_y_continuous(expand = c(0, 0)) + theme_classic() 
    # scale_fill_manual(values=  color_gradiant( c("blue","red"),1:length(levels(plot_data$taxa)) ) ) + 
    if ( manual_color_flag ) {
      p1 = p1 + scale_fill_manual( values = c("#000033","#3288bd","#66c2a5", "#fee08b", "#f46d43", "#d53e4f","#660000"),
                                   labels=new_legs, breaks = legs, name = 'Function' ) 
    } else {
      p1 = p1 + scale_fill_discrete( labels=new_legs, breaks = legs, name = 'Function')
    }
    
    #scale_fill_manual(values=c("#660000", "#d53e4f", "#f46d43","#fdae61", 
    #                           "#fee08b", "#ffffbf", "#e6f598","#abdda4")) #, 
    #                          "#66c2a5", "#3288bd", "#5e4fa2","#000033" )) 
    if (is.numeric( unlist( f_taxa[x_parm] ) ) )
      p1 = p1 + scale_x_continuous(expand = c(0, 0)) 
    else
    {
      if ( x_parm == 'Race' )
        p1 = p1 + scale_x_discrete(expand = c(0, 0), breaks=names, labels=new_names)
      else
        p1 = p1 + scale_x_discrete(expand = c(0, 0)) 
    }
  }
  
  return( list(p1,ttl, plot_data, f_taxa) )
}


# plots scatterplot af alpha dicersity from OTU file by filters and parameters.
plot_alpha_scatter = function(alpha_file, matadata_file, filter_parameters, met_cols, check_val, na_strings = c('no_data','_','NA','unknown'), rm_parameters = NA, rm_na=TRUE, one_color = F, print_pvals = F, sig_ast_flag = T )
{
  metadata = read.table(matadata_file,sep="\t", header=TRUE, na.strings = na_strings)
  alpha_div = read.table(alpha_file,sep="\t", header=TRUE, na.strings = na_strings)
  ## organizing metadata for alpha
  met = metadata[met_cols]
  alpha = data_add_metadata(alpha_div, met)
  
  # filtering data
  alpha = filter_data(alpha,met_parm)
  if ( length(rm_met_parm)!=0 )
    if ( !is.na(rm_parameters) )
      alpha = filter_data_rm(alpha,rm_parameters)
  
  ttl = sprintf('%s Alpha Diversity by %s',parameters_str( filter_parameters ),check_val)
  #alpha diversity scaterplot
  p1 = ggplot(alpha, aes(x=alpha[[check_val]], y=PD_whole_tree)) + geom_point(alpha = 0.6, size = 3, colour = "#194C99", shape = 19)  +
    stat_smooth(color = '#003366', size=1, fill="#6699CC", alpha = 0.3, method = loess, se=T) +
    theme(axis.text = element_text(colour = 'black'), axis.title=element_text(size=18), plot.title = element_text(size=22)) +
    xlab(check_val) + ylab('PD_whole_tree') + labs(title = ttl)  + theme_bw()# +
    #scale_shape_manual(values = 1)
  
  # p2 = ggplot2.boxplot(data=alpha, xName=check_val,yName='PD_whole_tree', groupName=check_val, groupColors=c('#999999','#E69F00') ) +  ##  for c_diff colors
  if (one_color)
    p2 = ggplot2.boxplot(data=alpha, xName=check_val,yName='PD_whole_tree', groupName='Sample_Location', groupColors = c('gold') ) # for one color
  else
    p2 = ggplot2.boxplot(data=alpha, xName=check_val,yName='PD_whole_tree', groupName=check_val ) 
  p2 = p2 + labs(title = ttl) + theme(axis.text.x = element_text( size=14), axis.text.y = element_text( size=14), axis.title = element_text(size = 16)) +
      xlab(check_val) +  ylab("Alpha Diversity") + theme_bw() + theme(legend.position="none") 
  
  if ( sig_ast_flag )
    p2 = add_significant_asterix_to_plot(p2, alpha[[check_val]], alpha[['PD_whole_tree']], print_pvals = print_pvals)
  return( list( p1, ttl, alpha, p2 ) )
}

# plots scatterplot af alpha dicersity from OTU file by filters and parameters.
plot_alpha_scatter_v2 = function(taxa_file, check_val, filter_parameters = NA, na_strings = c('no_data','_','NA','unknown'), rm_parameters = NA, rm_na=TRUE, file_flag = F, sig_ast_flag = T, jitter_flag = T, print_pvals = F, alpha_var = 'PD_whole_tree', y_lb = 'Alpha Diversity', fdr = 'bonferroni')
{
  if (file_flag)
    taxa = read.table(taxa_file,sep="\t", header=TRUE, na.strings = na_strings)
  else
    taxa = taxa_file
  taxa = taxa[ !is.na(taxa[[check_val]]), ]
  
  alpha = taxa
  # filtering data
  if ( length(rm_parameters)!=0 )
    if ( !is.na(rm_parameters) )
      alpha = filter_data_rm(alpha,rm_parameters)
  alpha = filter_data(alpha,filter_parameters, rmv_na = rm_na)
    
  ttl = sprintf('%s Alpha Diversity by %s',parameters_str( filter_parameters ),check_val)
  #alpha diversity scaterplot
  p1 = ggplot(alpha) + geom_point(aes(x=alpha[[check_val]], y=alpha[[alpha_var]]), alpha = 0.6, size = 3, colour = "#194C99", shape = 19)  +
    stat_smooth(aes(x=alpha[[check_val]], y=alpha[[alpha_var]]), color = '#003366', size=1, fill="#6699CC", alpha = 0.3, method = loess, se=T) +
    theme(axis.text = element_text(colour = 'black'), axis.title=element_text(size=18), plot.title = element_text(size=22)) +
    xlab(check_val) + ylab(y_lb) + labs(title = ttl)  + theme_bw()# +
  #scale_shape_manual(values = 1)
  
  # p2 = ggplot2.boxplot(data=alpha, xName=check_val,yName='PD_whole_tree', groupName=check_val, groupColors=c('#999999','#E69F00') ) +  ##  for c_diff colors
  p2 = ggplot(alpha) + geom_boxplot(aes(x=alpha[[check_val]],y = alpha[[alpha_var]], fill = alpha[[check_val]])) +
    xlab(check_val) +  ylab(y_lb)
  if ( jitter_flag )
    p2 = p2 + geom_jitter(aes(x=alpha[[check_val]], y=alpha[[alpha_var]]), alpha = 0.3) 
  p2 = format_fig(p2, x_angle = 45, set_ttl_flag = F) + theme(legend.position = 'none') + labs(title = ttl)
  
  # browser()
  if ( sig_ast_flag )
  {
    if ( fdr == 'bonferroni')
    {
      p2 = add_significant_asterix_to_plot_v2(p2, alpha[[check_val]], alpha[[alpha_var]], print_pvals = print_pvals)
    } else if ( fdr == 'BH')
    {
      p2 = add_significant_asterix_to_plot_BH(p2, alpha[[check_val]], alpha[[alpha_var]], print_pvals = print_pvals)
    }
  }
      return( list( p1, ttl, alpha, p2[[1]] ) )
}

# plots scatter plot if check_val is numeric and box plot if categoric. y should be numeric.
plot_variable_and_stats = function(taxa, check_val, y_val, sig_ast_flag = T, print_pvals = F, y_var = 'PD_whole_tree', y_lb = 'Alpha Diversity',boxplot_p_text_flag = T, boxplot_test_type = 'wilcox', jitter_flag = T, show_pval_as_asterix = T )
{
  if ( sum(!is.na(taxa[[check_val]])) == 0 )
    return( ggplot() )
  # scatter plot
  taxa = taxa[!is.na(taxa[[check_val]]),]
  if ( is.numeric( taxa[[ check_val ]] ) )
  {
    p = ggplot(taxa) + geom_point(aes_string(x=check_val, y=y_val)) +
      geom_smooth(aes_string(x=check_val, y=y_val), color = '#003366', size=1, fill="#6699CC", alpha = 0.3, method = loess, se=T) +
      xlab(check_val) + ylab(y_lb) + theme_bw() 
    if ( sig_ast_flag )
    {
      cor = cor.test(x=taxa[[check_val]], y=taxa[[y_val]], method = 'spearman')
      if (print_pvals)
      {
        library(grid)
        grob <- grobTree(textGrob(sprintf('spearman rho=%f, p=%f, n=%.0f', cor$estimate, cor$p.value, sum(!is.na(taxa[[check_val]]) )), x=0.1,  y=0.95, hjust=0) )
        p = p + annotation_custom(grob)
      }
    }
  }
  
  if (is.character( taxa[[ check_val ]] ) )
  {
    # adding n= to x labs
    temp = taxa[[check_val]]
    for ( i in 1:length(taxa[[check_val]] ) )
    {
      taxa[[check_val]][i] = sprintf('%s (n=%d)', temp[i], sum( temp == temp[i], na.rm = T ) )
      if (is.na(temp[i]))
        taxa[[check_val]][i] = sprintf('%s (n=%d)', temp[i], sum( is.na(temp) ) )
    }
    p = ggplot(taxa) + geom_boxplot(aes_string(x=check_val,y = y_val, fill = check_val), outlier.alpha = 0) +
      xlab(check_val) +  ylab(y_lb) +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none")
    if (jitter_flag)
      p = p + geom_jitter(aes_string(x=check_val, y=y_val), alpha = 0.3) 
    if ( sig_ast_flag )
    {
      if ( show_pval_as_asterix )
      {
        res = add_significant_asterix_to_plot_BH(p, as.factor(taxa[[check_val]]), taxa[[y_val]], print_pvals = print_pvals, test_type = boxplot_test_type, show_pval_as_asterix = show_pval_as_asterix)
      } else
      {
        res = add_significant_asterix_to_plot_BH(p, as.factor(taxa[[check_val]]), taxa[[y_val]], print_pvals = print_pvals, test_type = boxplot_test_type, show_pval_as_asterix = show_pval_as_asterix, label_size = 4, asterix_scale_var = 2.2)
      }
       p=res[[1]]
      t=res[[2]]
      if (dim(t)[1]>=1 & boxplot_p_text_flag == T)
      {
        t=t[t$p_val == min(t$p_val),]
        if (print_pvals)
        {
          library(grid)
          grob <- grobTree(textGrob(sprintf('%s min value %s vs %s:\np=%f, BH=%f', boxplot_test_type, t$var_1, t$var_2, t$p_val, t$BH), x=0.1,  y=0.95, hjust=0) )
          p = p + annotation_custom(grob)
        }
      }
    }
  }

  return( p )
}

# plots scatter plot if check_val is numeric and box plot if categoric. y should be numeric.
boxplot_variable_and_stats_permTest = function(taxa, check_val, y_val, sig_ast_flag = T, print_pvals = F, y_lb = 'Alpha Diversity',boxplot_p_text_flag = T, boxplot_test_type = 'perm', jitter_flag = T)
{
  if ( sum(!is.na(taxa[[check_val]])) == 0 | is.numeric( taxa[[ check_val ]] ))
    return( ggplot() )
  # scatter plot
  taxa = taxa[!is.na(taxa[[check_val]]),]

  # if (is.character( taxa[[ check_val ]] ) )
  {
    # adding n= to x labs
    temp = taxa[[check_val]]
    for ( i in 1:length(taxa[[check_val]] ) )
    {
      taxa[[check_val]][i] = sprintf('%s (n=%d)', temp[i], sum( temp == temp[i], na.rm = T ) )
      if (is.na(temp[i]))
        taxa[[check_val]][i] = sprintf('%s (n=%d)', temp[i], sum( is.na(temp) ) )
    }
    p = ggplot(taxa) + geom_boxplot(aes_string(x=check_val,y = y_val, fill = check_val), outlier.alpha = 0) +
      xlab(check_val) +  ylab(y_lb) + 
      theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none")
    if (jitter_flag)
      p = p+geom_jitter(aes_string(x=check_val, y=y_val), alpha = 0.3) 
    if ( sig_ast_flag )
    {
      res = add_significant_asterix_to_plot_BH_v2(p, as.factor(taxa[[check_val]]), taxa[[y_val]], print_pvals = print_pvals, test_type = boxplot_test_type)
      p=res[[1]]
      t=res[[2]]
      if (dim(t)[1]>=1 & boxplot_p_text_flag == T)
      {
        t=t[t$p_val == min(t$p_val),]
        if (print_pvals)
        {
          library(grid)
          grob <- grobTree(textGrob(sprintf('%s min value %s vs %s:\np=%f, BH=%f', boxplot_test_type, t$var_1, t$var_2, t$p_val, t$BH), x=0.1,  y=0.95, hjust=0) )
          p = p + annotation_custom(grob)
        }
      }
    }
  }
  
  return( p )
}


# removing taxa with name containing a dot (repeats of same genus with taxonomic differnce)
remove_dot_taxa = function(FC_vlas, taxa_pos = 'taxa')
{
  pos = vector(mode='logical',length = length(FC_vlas[[taxa_pos]]) )
  for (i in  1:length(FC_vlas[[taxa_pos]]))
  {
    pos[i] = !grepl('[.]',FC_vlas[[taxa_pos]][i]) #tells if there is a dit in the name
  }
  FC_vlas = FC_vlas[pos,]
  FC_vlas = droplevels(FC_vlas)
}

# makes color gradient of between wanted 2 colors
# values = 1:7 for example for a 7-color gradient
color_gradiant <- function(colors, values) 
{
   v <- (values - min(values))/diff(range(values))
   x <- colorRamp(colors)(v)
   rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}

# formats a figures text size, title, axis angle etc.
format_fig = function(fig, lab = NA, title_size = 20, text_size = 16, axis_size = 12, x_angle = NA, set_ttl_flag = T, bw_flag = T)
{
  fig = fig + 
    theme(text = element_text(size=text_size),
          axis.text.x = element_text(size=axis_size))
  if ( bw_flag )
    fig = fig + theme_bw()
  if ( set_ttl_flag )
    fig = fig + labs(title = lab) + theme(plot.title = element_text(hjust=0, size = title_size, vjust = 1))
  if ( !is.na(x_angle) )
  {
    fig = fig + theme( axis.text.x = element_text(angle = x_angle, vjust = 1, hjust = 1))
    if ( x_angle ==90 )
      fig = fig + theme( axis.text.x = element_text(angle = x_angle, vjust = 0.5, hjust = 1))
  }
  return(fig)
}

# finds beta diversity distance between two sample. uses names of two samples and a beta diversity table
two_samps_beta_dist = function(smp1, smp2, bDiv)
{
  ids1 = as.character(bDiv[,1])
  ids2 = as.character(bDiv$V1)
  
  pos1 =  which(smp1 == ids1)
  pos2 =  which(smp2 == ids2)
  
  return(as.numeric(as.character(bDiv[pos2, pos1])))
}

# finds positions of strings in to_search (character vector) inside of to_search_in (a bigger character vector)
strings_in_strings_pos = function(to_search, to_search_in)
{
  pos = vector( mode = 'numeric',length = length(to_search) )
  for (i in 1:length(to_search))
  {
    if ( to_search[i] %in% to_search_in )
      pos[i] = which(to_search[i] == to_search_in)
    # else
      # print(to_search[i])
  }
  return(pos)
}

# finds beta diversity distances values between two sample lists. uses names of the samples and a beta diversity table
samps_beta_dist = function(smps1, smps2, bDiv)
{
  smps1 = as.character(smps1)
  smps2 = as.character(smps2)
  
  ids1 = as.character(bDiv[,1])
  ids2 = as.character(bDiv$V1)
  
  pos1 = strings_in_strings_pos(smps1, ids1)
  pos2 = strings_in_strings_pos(smps2, ids2)
  
  x=as.numeric(c(t(bDiv[pos2, pos1])))
  x=x[x!=0]
  
  return(x)
}

# makes a boxplot of beta diversity distance whithin groups by check_val grouping.
# return figure and dataframe of distance values for groups.
# get beta diversity file, taxa dataframe, name of parameter to chack.
# also (optional): list of filtering parameters, weight (onlt matters for title) and na string for file read 
claculate_bDiv_boxplot = function(b_file, taxa, check_val, flt_prm = NA, weight = 'unweighted', na_str = c('no_data','_','NA','unknown', 'other'), check_between_flag = FALSE, sig_ast_flag = T,sort_by_Cohort_flag = F, test_type = 't_test', ord = c(5,4,3,6,7,1,2) )
{
  # read beta div file
  bDiv = read.table(b_file,sep="\t", header=FALSE, na.strings = na_str)
  # filter data by parameters
  taxa =  filter_data(taxa,flt_prm)
  # get title by parameters
  flt_ttl = parameters_str(flt_prm)
  
  lvls = levels(taxa[[check_val]])
  if ( is.numeric(taxa[[check_val]]) )
    lvls = as.character(unique(taxa[[check_val]]))
    
  value = c()
  variable = c()
  if (check_between_flag) 
  {
    for ( i in 1:length(lvls) ) 
    {
      for (j in 1:i)
      {
        lvl_1 = lvls[i]
        lvl_2 = lvls[j]
        
        samps1 = taxa$SampleID[taxa[[check_val]]==lvl_1]
        samps2 = taxa$SampleID[taxa[[check_val]]==lvl_2]
        
        beta_vals = samps_beta_dist(samps1, samps2, bDiv)
        value = c(value, beta_vals)
        
        lab = sprintf('%s_v_%s',lvl_1, lvl_2)
        variable = c(variable, rep(lab, length(beta_vals)))
      }
    }
  } else 
  {
    for (lvl in lvls)
    {
      samps = taxa$SampleID[taxa[[check_val]]==lvl]
      
      beta_vals = samps_beta_dist(samps, samps, bDiv)
      value = c(value, beta_vals)
      variable = c(variable, rep(lvl, length(beta_vals)))
    }
  }
  
  data = data.frame(value = value, variable = variable)
  
  ttl = sprintf('Beta Diversity %s Distance %s', check_val, flt_ttl)
  if ( check_between_flag )
    ttl = sprintf('%s between',ttl)
  
  if ( sort_by_Cohort_flag )
  {
    a=levels(data$variable)
    data$variable <- factor(data$variable, levels = a[ord] )
  }
  
  p = ggplot(data) + geom_boxplot(aes(x= factor(variable), y=value, fill = factor(data$variable))) + theme_bw() + 
    xlab(check_val) +  ylab(sprintf("%s UniFrac Distance",weight))  + labs(title = ttl) +
    theme(legend.position = 'none') + theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  if ( sig_ast_flag )
    p = add_significant_asterix_to_plot(p, factor(data$variable), data$value, print_pvals = T, test_type = test_type)
  
  return(list(p,data, ttl))
}

# makes a boxplot of beta diversity distance whithin groups by check_val grouping.
# return figure and dataframe of distance values for groups.
# get beta diversity file, taxa dataframe, name of parameter to chack.
# also (optional): list of filtering parameters, weight (onlt matters for title) and na string for file read 
# v2 = no repeats of samp1-samp2, samp2-samp1 etc. in within group comparision
claculate_bDiv_boxplot_v2 = function(b_file, taxa, check_val, flt_prm = NA, weight = 'unweighted', na_str = c('no_data','_','NA','unknown', 'other'), check_between_flag = FALSE, sig_ast_flag = T,sort_by_Cohort_flag = F, test_type = 't_test', ord = c(5,4,3,6,7,1,2) )
{
  # read beta div file
  bDiv = read.table(b_file,sep="\t", header=FALSE, na.strings = na_str)
  # filter data by parameters
  taxa =  filter_data(taxa,flt_prm)
  # get title by parameters
  flt_ttl = parameters_str(flt_prm)
  
  lvls = levels(taxa[[check_val]])
  if ( is.numeric(taxa[[check_val]]) )
    lvls = as.character(unique(taxa[[check_val]]))
  
  value = c()
  variable = c()
  if (check_between_flag) 
  {
    for ( i in 1:length(lvls) ) 
    {
      for (j in 1:i)
      {
        lvl_1 = lvls[i]
        lvl_2 = lvls[j]
        
        samps1 = taxa$SampleID[taxa[[check_val]]==lvl_1]
        samps2 = taxa$SampleID[taxa[[check_val]]==lvl_2]
        
        beta_vals = samps_beta_dist(samps1, samps2, bDiv)
        value = c(value, beta_vals)
        
        lab = sprintf('%s_v_%s',lvl_1, lvl_2)
        variable = c(variable, rep(lab, length(beta_vals)))
      }
    }
  } else 
  {
    for (lvl in lvls)
    {
      # samps = taxa$SampleID[taxa[[check_val]]==lvl]
      pos = which( taxa[[check_val]]==lvl )
      
      beta_vals = c()
      for ( i in 1:length(pos) )
      {
        for ( j in 1:i )
        {
          if (i!=j)
            beta_vals = c(beta_vals, bDiv[i,j])
        }
      }
      # beta_vals = samps_beta_dist(samps, samps, bDiv)
      
      value = c(value, beta_vals)
      variable = c(variable, rep(lvl, length(beta_vals)))
    }
  }
  
  data = data.frame(value = value, variable = variable)
  
  ttl = sprintf('Beta Diversity %s Distance %s', check_val, flt_ttl)
  if ( check_between_flag )
    ttl = sprintf('%s between',ttl)
  
  if ( sort_by_Cohort_flag )
  {
    a=levels(data$variable)
    data$variable <- factor(data$variable, levels = a[ord] )
  }
  
  p = ggplot(data) + geom_boxplot(aes(x= factor(variable), y=value, fill = factor(data$variable))) + theme_bw() + 
    xlab(check_val) +  ylab(sprintf("%s UniFrac Distance",weight))  + labs(title = ttl) +
    theme(legend.position = 'none') + theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  if ( sig_ast_flag )
    p = add_significant_asterix_to_plot(p, factor(data$variable), data$value, print_pvals = T, test_type = test_type)
  
  return(list(p,data, ttl))
}

# makes a boxplot of beta diversity distance whithin groups by check_val grouping.
# return figure and dataframe of distance values for groups.
# get beta diversity file, taxa dataframe, name of parameter to chack.
# also (optional): list of filtering parameters, weight (onlt matters for title) and na string for file read 
claculate_bDiv_boxplot_allBetween = function(b_file, taxa, check_val, flt_prm = NA, weight = 'unweighted', na_str = c('no_data','_','NA','unknown', 'other'), sig_ast_flag = T,sort_by_Cohort_flag = F, test_type = 't_test', ord = c(5,4,3,6,7,1,2) )
{
  # read beta div file
  bDiv = read.table(b_file,sep="\t", header=FALSE, na.strings = na_str)
  # filter data by parameters
  taxa =  filter_data(taxa,flt_prm)
  # get title by parameters
  flt_ttl = parameters_str(flt_prm)
  
  lvls = levels(taxa[[check_val]])
  if ( is.numeric(taxa[[check_val]]) )
    lvls = as.character(unique(taxa[[check_val]]))
  
  value = c()
  variable = c()
  if (check_between_flag) 
  {
    for ( i in 1:length(lvls) ) 
    {
      for (j in 1:i)
      {
        lvl_1 = lvls[i]
        lvl_2 = lvls[j]
        
        samps1 = taxa$SampleID[taxa[[check_val]]==lvl_1]
        samps2 = taxa$SampleID[taxa[[check_val]]==lvl_2]
        
        beta_vals = samps_beta_dist(samps1, samps2, bDiv)
        value = c(value, beta_vals)
        
        lab = sprintf('%s_v_%s',lvl_1, lvl_2)
        variable = c(variable, rep(lab, length(beta_vals)))
      }
    }
  } else 
  {
    for (lvl in lvls)
    {
      samps = taxa$SampleID[taxa[[check_val]]==lvl]
      
      beta_vals = samps_beta_dist(samps, samps, bDiv)
      value = c(value, beta_vals)
      variable = c(variable, rep(lvl, length(beta_vals)))
    }
  }
  
  data = data.frame(value = value, variable = variable)
  
  ttl = sprintf('Beta Diversity %s Distance %s', check_val, flt_ttl)
  if ( check_between_flag )
    ttl = sprintf('%s between',ttl)
  
  if ( sort_by_Cohort_flag )
  {
    a=levels(data$variable)
    data$variable <- factor(data$variable, levels = a[ord] )
  }
  
  p = ggplot(data) + geom_boxplot(aes(x= factor(variable), y=value, fill = factor(data$variable))) + theme_bw() + 
    xlab(check_val) +  ylab(sprintf("%s UniFrac Distance",weight))  + labs(title = ttl) +
    theme(legend.position = 'none') + theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  if ( sig_ast_flag )
    p = add_significant_asterix_to_plot(p, factor(data$variable), data$value, print_pvals = T, test_type = test_type)
  
  return(list(p,data, ttl))
}


# calculating microbial dysbiosis index, as in RISK and longitues studies.
calc_MD_index = function(taxa, positive_taxa = 
                           c('Enterobacteriaceae','Pasteurellaceae','Fusobacteriaceae','Neisseriaceae','Veillonellaceae','Gemellaceae'), 
                         negative_taxa = c('Bacteroidales','Clostridiales','Erysipelotrichaceae','Bifidobacteriaceae'), 
                         bad_negative = 'Veillonellaceae', check_level = 6, type = 'MD')
{
  # get taxa names for positive/negative taxa position use later
  tx_names = names(taxa)
  
  # positive and negative taxa search string
  positive_tx = paste(positive_taxa,sep='',collapse ='|')
  negative_tx = paste(negative_taxa,sep='',collapse ='|')
  
  # position of positive and negative taxa in taxa dataframe
  positive_pos = grepl(positive_tx,tx_names) & check_id_level(tx_names) == check_level
  negative_pos = intersect(grep(negative_tx,tx_names), grep(bad_negative,tx_names, invert = T) )
  negative_pos = intersect(negative_pos, which(check_id_level(tx_names) == check_level) )
  
  # calculate microbial dysbiosis index
  l = dim(taxa)[1]
  dys_inx = vector( mode = 'numeric', length = l)
  for ( i in 1:l )
  {
    sample_pos = i
    if ( type == 'MD' ) {
      dys_inx[i] = log10( sum(taxa[sample_pos,positive_pos])/sum(taxa[sample_pos,negative_pos]) )
    } else if ( type == 'inc' ) {
      dys_inx[i] = sum( taxa[sample_pos,positive_pos] )
    } else if ( type == 'dec' ) {
      dys_inx[i] = sum( taxa[sample_pos,negative_pos] ) 
    } else
      print('type not one of recognized types (MD, inc or dec)')
  }
  dys_inx[is.infinite(dys_inx)]=NA
  
  return(dys_inx)
}

# plots dysbiosis index boxplot by check_val, with significant asterix if wanted
plot_dysbiosis_index = function(taxa, check_val, x_angle = 45, sig_ast_flag = T, print_pvals = F, test_type = 't_test')
{
  taxa = taxa[ !is.na(taxa[[check_val]]), ]
  
  ttl = sprintf('Dysbiosis Index Boxplot by %s', check_val)
  
  p = ggplot(taxa) + geom_boxplot(aes(x=taxa[[check_val]], y=dysbiosis_index, fill = taxa[[check_val]])) + 
    xlab(check_val) # + geom_jitter(aes(x=taxa[[check_val]], y=dysbiosis_index), alpha = 0.3) 
  p = format_fig(p, x_angle = 45, set_ttl_flag = F) + theme(legend.position = 'none') + labs(title = ttl)
  
  if ( sig_ast_flag )
  {
    var = taxa[[check_val]]
    val = taxa$dysbiosis_index
    p = add_significant_asterix_to_plot(p, var, val, print_pvals = print_pvals, test_type = test_type)
  }
  
  return(p)
}

# change a numeric variable to bimary Yes/No by cutoff
numeric_to_binary = function(var, cutoff, low_str = 'Y', high_str = 'N')
{
  new_var = vector(mode = 'character', length = length(var) )
  for ( i in 1:length(var) )
  {
    if ( is.na(var[i] ))
    {
      new_var[i] = NA
    } else
      new_var[i] = if (var[i]<cutoff) low_str else high_str
  }
  return(new_var)
}

# change a numeric variable to discrete 3 options variable
numeric_to_tertiary = function(var, cutoff, high_cutoff, low_str = 'Y', med_str = 'M',high_str = 'N')
{
  new_var = vector(mode = 'character', length = length(var) )
  for ( i in 1:length(var) )
  {
    if ( is.na(var[i] ))
    {
      new_var[i] = NA
    } else
      new_var[i] = if (var[i]<cutoff) low_str else if (var[i]<high_cutoff) med_str else high_str
  }
  return(new_var)
}


# find fitting asterix to represent the p_value (usually in plots)
find_pVal_asterix = function(p_val)
{
  if (p_val <= 0.001)
    { return('***') }
  if (p_val <= 0.01)
    { return('**') }
  if (p_val <= 0.05)
    { return('*') }
  return('ns')
}

# adds significant asterix by p. value to plot, by Bonferroni correncted t/wilcox test 
# p = the plot to add on. has to be built as ggplot(data = ?), with aes in geom_boxplot etc.
# var = the variable to divide by. needs to be factor
# val = the values
# p_val_cutoff = cuoff of p.value to show differnce form. default is 0.05
# line_width = width of the line. astetic thing
add_significant_asterix_to_plot = function(p, var, val, p_val_cutoff = 0.05, line_width = 1, print_pvals = F, log_scale_flag=F, test_type = 't_test')
{
  all_vars = levels(var)
  
  # calculate the number of comparisons for multiple comparisons Bonferroni correction
  l = length( all_vars )
  fdr_fix = (l*(l-1))/2
  
  scale = sd(val, na.rm = T)/6
  # if ( log_scale_flag )
  #   scale = 2^scale
  
  max_p = -1
  max_pos = max(val, na.rm = T) + scale
  # going over all comparisons between differnt variable (no repeats)
  for ( i in 1:l )
  {
    if ( i-1 != 0 ) # because 1:0 is 1 0 and not nothing in R
    {
      for ( j in 1:(i-1) )
      {
        var_1 = all_vars[i]
        var_2 = all_vars[j]
        
        pos1 = var == var_1 & !is.na(var)
        pos2 = var == var_2 & !is.na(var)
        
        # if there is enough data to check
        res = list()
        if ( sum(pos1)>1 & sum(pos2)>1 )
        {
          # calculate p.vlaue
          if ( test_type == 't_test') {
            res = t.test(x = val[pos1], y = val[pos2], na.rm = T)
          } else if ( test_type == 'wilcox' ) {
            res = wilcox.test(x = val[pos1], y = val[pos2], na.rm = T)
          } else {print ('test type does not fit to "t_test" or "wilcox')}
          if (print_pvals)
          {
            print(sprintf('%s vs %s p.val = %f', var_1, var_2, res$p.value * fdr_fix ))
            if ( res$p.value* fdr_fix < p_val_cutoff )
              if ( res$p.value* fdr_fix > max_p )
                max_p = res$p.value* fdr_fix
          }
        } else {res$p.value = 1}
        # res = wilcox.test(x = val[var == var_1], y = val[var == var_2])
        p_val = res$p.value * fdr_fix
        if ( p_val < p_val_cutoff ) # if p.value is significant enough, add asterix
        {
          # paramters
          astx_pos = max_pos
          start_pos = i
          end_pos = j
          lbl = find_pVal_asterix(p_val)
          
          # building neccecary data frame
          x_pos = c(start_pos, start_pos,end_pos,end_pos)
          y_pos = c(astx_pos, astx_pos+scale*1, astx_pos+scale*1, astx_pos)
          df2 <- data.frame(x_pos = x_pos, y_pos = y_pos )
          
          # adding the asterix to plot
          p = p + geom_path(data = df2, aes(x = x_pos, y = y_pos),size = line_width) + 
            annotate("text", x = mean( c(start_pos, end_pos) ), y = astx_pos+scale*1.3, label = lbl, size = 8) 
          
          # updating maximal position on plot for next time
          max_pos = max_pos + scale*2.5
          if (log_scale_flag)
            max_pos = max_pos + scale*max_pos
        }
        
      }
    }
  }
  if (print_pvals)
    print(sprintf('maximal significant p.val = %e', max_p ))
  return(p)
}

# adds significant asterix by p. value to plot, by Bonferroni correncted t/wilcox test 
# p = the plot to add on. has to be built as ggplot(data = ?), with aes in geom_boxplot etc.
# var = the variable to divide by. needs to be factor
# val = the values
# p_val_cutoff = cuoff of p.value to show differnce form. default is 0.05
# line_width = width of the line. astetic thing
# v2 = return list with plot and maximal significant p value
add_significant_asterix_to_plot_v2 = function(p, var, val, p_val_cutoff = 0.05, line_width = 1, print_pvals = F, log_scale_flag=F, test_type = 't_test',label_size = 8, scale_var = 6, asterix_scale_var = 1.3)
{
  all_vars = levels(var)
  
  # calculate the number of comparisons for multiple comparisons Bonferroni correction
  l = length( all_vars )
  fdr_fix = (l*(l-1))/2
  
  scale = sd(val, na.rm = T)/scale_var
  # if ( log_scale_flag )
  #   scale = 2^scale
  
  max_p = -1
  max_pos = max(val, na.rm = T) + scale
  # going over all comparisons between differnt variable (no repeats)
  for ( i in 1:l )
  {
    if ( i-1 != 0 ) # because 1:0 is 1 0 and not nothing in R
    {
      for ( j in 1:(i-1) )
      {
        var_1 = all_vars[i]
        var_2 = all_vars[j]
        
        pos1 = var == var_1 & !is.na(var)
        pos2 = var == var_2 & !is.na(var)
        
        # if there is enough data to check
        res = list()
        if ( sum(pos1)>1 & sum(pos2)>1 )
        {
          # calculate p.vlaue
          if ( test_type == 't_test') {
            ##                            big problem!!!!!!!!!!!
            # tryCatch( res = t.test(x = val[pos1], y = val[pos2], na.rm = T),
            #           error = function(e){print(sprintf('t_test error %s',e$message) )}, 
            #           finally = {res$p.value = 1})
            res = t.test(x = val[pos1], y = val[pos2], na.rm = T)
          } else if ( test_type == 'wilcox' ) {
            res = wilcox.test(x = val[pos1], y = val[pos2], na.rm = T)
          } else {print ('test type does not fit to "t_test" or "wilcox')}
          if (print_pvals)
          {
            print(sprintf('%s vs %s p.val = %f * %i = %f', var_1, var_2, res$p.value, fdr_fix, res$p.value * fdr_fix ))
          }
          if ( !is.na(res$p.value) & res$p.value* fdr_fix < p_val_cutoff )
            if ( res$p.value* fdr_fix > max_p )
              max_p = res$p.value* fdr_fix
        } else {res$p.value = 1}
        # res = wilcox.test(x = val[var == var_1], y = val[var == var_2])
        p_val = res$p.value * fdr_fix
        if ( !is.na(p_val))
          if ( p_val < p_val_cutoff ) # if p.value is significant enough, add asterix
          {
            # paramters
            astx_pos = max_pos
            start_pos = i
            end_pos = j
            lbl = find_pVal_asterix(p_val)
            
            # building neccecary data frame
            x_pos = c(start_pos, start_pos,end_pos,end_pos)
            y_pos = c(astx_pos, astx_pos+scale*1, astx_pos+scale*1, astx_pos)
            df2 <- data.frame(x_pos = x_pos, y_pos = y_pos )
            
            # adding the asterix to plot
            p = p + geom_path(data = df2, aes(x = x_pos, y = y_pos),size = line_width) + 
              annotate("text", x = mean( c(start_pos, end_pos) ), y = astx_pos+scale*asterix_scale_var, label = lbl, size = label_size) 
            
            # updating maximal position on plot for next time
            max_pos = max_pos + scale*2.5
            if (log_scale_flag)
              max_pos = max_pos + scale*max_pos
          }
        
      }
    }
  }
  if (print_pvals)
    print(sprintf('maximal significant p.val = %e', max_p ))
  if (max_p == -1)
  {
    max_p = 1
  }
  return( list(p,max_p) )
}



# adds significant asterix by p. value to plot, 
# by Benjamini–Hochberg procedure (not Bonferroni correncted here) and wilcox test 
# p = the plot to add on. has to be built as ggplot(data = ?), with aes in geom_boxplot etc.
# var = the variable to divide by. needs to be factor
# val = the values
# p_val_cutoff = cuoff of p.value to show differnce form. default is 0.05
# line_width = width of the line. astetic thing
add_significant_asterix_to_plot_BH = function(p, var, val, p_val_cutoff = 0.05, line_width = 1, print_pvals = F, log_scale_flag=F, test_type = 't_test',label_size = 8, scale_var = 6, asterix_scale_var = 1.3, manual_label = '', show_pval_as_asterix = T)
{
  all_vars = levels(var)
  
  # calculate the number of comparisons
  l = length( all_vars )
  
  scale = sd(val, na.rm = T)/scale_var
  # if ( log_scale_flag )
  #   scale = 2^scale
  
  max_pos = max(val, na.rm = T) + scale
  # going over all comparisons between differnt variable (no repeats)
  # this repeat is to calculate p vals
  
  p_vals = c()
  vars_1 = c()
  vars_2 = c()
  is = c()
  js = c()
  for ( i in 1:l )
  {
    if ( i-1 != 0 ) # because 1:0 is 1 0 and not nothing in R
    {
      for ( j in 1:(i-1) )
      {
        var_1 = all_vars[i]
        var_2 = all_vars[j]
        
        pos1 = var == var_1 & !is.na(var)
        pos2 = var == var_2 & !is.na(var)
        
        # if there is enough data to check
        res = list()
        if ( sum(pos1)>1 & sum(pos2)>1 )
        {
          # calculate p.vlaue
          if ( test_type == 't_test') {
            ##                            big problem!!!!!!!!!!!
            # tryCatch( res = t.test(x = val[pos1], y = val[pos2], na.rm = T),
            #           error = function(e){print(sprintf('t_test error %s',e$message) )}, 
            #           finally = {res$p.value = 1})
            res = t.test(x = val[pos1], y = val[pos2], na.rm = T)
          } else if ( test_type == 'wilcox' ) {
            res = wilcox.test(x = val[pos1], y = val[pos2], na.rm = T)
          } else {print ('test type does not fit to "t_test" or "wilcox')}
        } else {res$p.value = NA}

        p_val = res$p.value 
        if ( !is.na(p_val) ) # if a legit p value (enough data to make a comparison, etc.), update it
        {
          p_vals = c(p_vals, p_val)
          vars_1 = c(vars_1, var_1)
          vars_2 = c(vars_2, var_2)
          is = c(is, i)
          js = c(js, j)
        }
      }
    }
  }
  # calculating Benjamini–Hochberg
  BH = p.adjust(p_vals,method="BH") 
  
  # going over all comparisons between differnt variable (no repeats)
  for ( i in 1:l )
  {
    if ( i-1 != 0 ) # because 1:0 is 1 0 and not nothing in R
    {
      for ( j in 1:(i-1) )
      {
        var_1 = all_vars[i]
        var_2 = all_vars[j]
        
        pos = which( is == i & js == j )
        
        if ( length(pos) >1 ) { print('ERROR! check the p value assinged'); return() }
        if ( length(pos) != 0 ) # if there is enough data to check
        {
          if (print_pvals)
          {
            print(sprintf('%s vs %s p.val = %f , BH = %f', var_1, var_2, p_vals[pos], BH[pos] ))
          }
          if ( BH[pos] < p_val_cutoff ) # if p.value is significant enough, add asterix
          {
            # paramters
            astx_pos = max_pos
            start_pos = i
            end_pos = j
            if (manual_label == '')
            {
              lbl = find_pVal_asterix(BH[pos])
              if (!show_pval_as_asterix)
              {
                lbl = ifelse( BH[pos] < 0.009, sprintf('Q = %.2e',BH[pos]), sprintf('Q = %.4s',BH[pos]) )
              }
            }  else 
              lbl = manual_label
            
            # building neccecary data frame
            x_pos = c(start_pos, start_pos,end_pos,end_pos)
            y_pos = c(astx_pos, astx_pos+scale*1, astx_pos+scale*1, astx_pos)
            df2 <- data.frame(x_pos = x_pos, y_pos = y_pos )
            
            # adding the asterix to plot
            p = p + geom_path(data = df2, aes(x = x_pos, y = y_pos),size = line_width) + 
              annotate("text", x = mean( c(start_pos, end_pos) ), y = astx_pos+scale*asterix_scale_var, label = lbl, size = label_size) 
            
            # updating maximal position on plot for next time
            max_pos = max_pos + scale*2.5
            if (log_scale_flag)
              max_pos = max_pos + scale*max_pos
          }
        }
      }
    }
  }
  # creating a dataframe with the p value the BH corrected values per varibles 
  pvals_df = data.frame( var_1 = vars_1, var_2 = vars_2, p_val = p_vals, BH = BH )
  return( list(p, pvals_df) )
}

# get a dataframe fit to add significance datat to plots using ggsignif package
# @test_type : possible test type are: t_test, wilcox, perm
# @fdr_flag : fdr using BH
get_signif_data_for_plot = function(var, val, test_type = 'wilcox', 
                                    fdr_flag = T, p_val_cutoff = 0.05,  
                                    show_pval_as_asterix = T,
                                    print_pvals = F, log_scale_flag=F, 
                                    scale_var = 6, asterix_scale_var = 1.3 )
{
  all_vars = unique(var)

  # calculate the number of comparisons
  l = length( all_vars )
  
  scale = sd(val, na.rm = T)/scale_var
  max_pos = max(val, na.rm = T) + scale
  
  # going over all comparisons between different variable (no repeats)
  # this repeat is to calculate p values
  
  p_vals = c()
  vars_1 = c()
  vars_2 = c()
  is = c()
  js = c()
  for ( i in 1:l )
  {
    if ( i-1 != 0 ) # because 1:0 is 1 0 and not nothing in R
    {
      for ( j in 1:(i-1) )
      {
        var_1 = all_vars[i]
        var_2 = all_vars[j]
        
        pos1 = var == var_1 & !is.na(var)
        pos2 = var == var_2 & !is.na(var)
        
        # if there is enough data to check
        res = list()
        if ( sum(pos1)>1 & sum(pos2)>1 )
        {
          # calculate p.vlaue
          if ( test_type == 't_test') {
            res = t.test(x = val[pos1], y = val[pos2], na.rm = T)
          } else if ( test_type == 'wilcox' ) {
            res = wilcox.test(x = val[pos1], y = val[pos2], na.rm = T)
          } else if ( test_type == 'perm' )
          {
            res = perm_test(x = val[pos1], y = val[pos2])
          } else {print ('test type does not fit to "t_test", "wilcox" or "perm". ')}
        } else {res$p.value = NA}
        
        p_val = res$p.value 
        if ( !is.na(p_val) ) # if a legit p value (enough data to make a comparison, etc.), update it
        {
          p_vals = c(p_vals, p_val)
          vars_1 = c(vars_1, var_1)
          vars_2 = c(vars_2, var_2)
          is = c(is, i)
          js = c(js, j)
        }
      }
    }
  }
  # calculating Benjamini–Hochberg
  BH = p.adjust(p_vals,method="BH") 
  
  res_df = data.frame(pval = p_vals, fdr = BH, var1 = vars_1, var2 = vars_2,
                      xmin =is, xmax=js, annotation = 0, y_position = 0, res_value = 2)
  # browser()
  for (i in 1:dim(res_df)[1] )
  {
    if (print_pvals)
    {
      print(sprintf('%s vs %s p.val = %f , BH = %f', res_df$var1[i], res_df$var2[i], 
                    res_df$pval[i], res_df$fdr[i] ))
    }
    res_df$res_value[i] = ifelse(fdr_flag, res_df$fdr[i], res_df$pval[i])
    if ( res_df$res_value[i] <= p_val_cutoff ) # if p.value is significant enough, add asterix
    {
      lbl = res_df$res_value[i]
      if ( show_pval_as_asterix )
      {
        lbl = find_pVal_asterix(lbl)
      } else
      {
        lbl = ifelse( lbl < 0.009, sprintf('Q = %.2e',lbl), sprintf('Q = %.4s',lbl) )
        if ( l == 2 | !fdr_flag )  
          lbl = gsub('Q','P',lbl)
      }
      res_df$annotation[i] = lbl
      
      res_df$y_position[i] = max_pos
      # updating maximal position on plot for next time
      max_pos = max_pos + scale*2.5
      if (log_scale_flag)
        max_pos = max_pos + scale*max_pos
    }
  }
  res_df = res_df[res_df$res_value <= p_val_cutoff,]
  return( res_df )
}

# adds significant asterix by p. value to plot, 
# by Benjamini–Hochberg procedure (not Bonferroni correncted here) and wilcox test 
# p = the plot to add on. has to be built as ggplot(data = ?), with aes in geom_boxplot etc.
# var = the variable to divide by. needs to be factor
# val = the values
# p_val_cutoff = cuoff of p.value to show differnce form. default is 0.05
# line_width = width of the line. astetic thing
# v2 - added new test type - self made permutation (ranked mean difference, 1000+ permutaion, x>= my value ,x+1/x+1) test, as Amnon recommended.
add_significant_asterix_to_plot_BH_v2 = function(p, var, val, p_val_cutoff = 0.05, line_width = 1, print_pvals = F, log_scale_flag=F, test_type = 'perm',label_size = 8, scale_var = 6, asterix_scale_var = 1.3, manual_label = '', show_pval_as_asterix = T)
{
  all_vars = levels(var)
  
  # calculate the number of comparisons
  l = length( all_vars )
  
  scale = sd(val, na.rm = T)/scale_var
  # if ( log_scale_flag )
  #   scale = 2^scale
  
  max_pos = max(val, na.rm = T) + scale
  # going over all comparisons between differnt variable (no repeats)
  # this repeat is to calculate p vals
  
  p_vals = c()
  vars_1 = c()
  vars_2 = c()
  is = c()
  js = c()
  for ( i in 1:l )
  {
    if ( i-1 != 0 ) # because 1:0 is 1 0 and not nothing in R
    {
      for ( j in 1:(i-1) )
      {
        var_1 = all_vars[i]
        var_2 = all_vars[j]
        
        pos1 = var == var_1 & !is.na(var)
        pos2 = var == var_2 & !is.na(var)
        
        # if there is enough data to check
        res = list()
        if ( sum(pos1)>1 & sum(pos2)>1 )
        {
          # calculate p.vlaue
          if ( test_type == 't_test') {
            ##                            big problem!!!!!!!!!!!
            # tryCatch( res = t.test(x = val[pos1], y = val[pos2], na.rm = T),
            #           error = function(e){print(sprintf('t_test error %s',e$message) )}, 
            #           finally = {res$p.value = 1})
            res = t.test(x = val[pos1], y = val[pos2], na.rm = T)
          } else if ( test_type == 'wilcox' ) {
            res = wilcox.test(x = val[pos1], y = val[pos2], na.rm = T)
          } else if ( test_type == 'perm' )
          {
            res = perm_test(x = val[pos1], y = val[pos2])
          } else {print ('test type does not fit to "t_test", "wilcox" or "perm". ')}
        } else {res$p.value = NA}
        
        p_val = res$p.value 
        if ( !is.na(p_val) ) # if a legit p value (enough data to make a comparison, etc.), update it
        {
          p_vals = c(p_vals, p_val)
          vars_1 = c(vars_1, var_1)
          vars_2 = c(vars_2, var_2)
          is = c(is, i)
          js = c(js, j)
        }
      }
    }
  }
  # calculating Benjamini–Hochberg
  BH = p.adjust(p_vals,method="BH") 
  
  # going over all comparisons between differnt variable (no repeats)
  for ( i in 1:l )
  {
    if ( i-1 != 0 ) # because 1:0 is 1 0 and not nothing in R
    {
      for ( j in 1:(i-1) )
      {
        var_1 = all_vars[i]
        var_2 = all_vars[j]
        
        pos = which( is == i & js == j )
        
        if ( length(pos) >1 ) { print('ERROR! check the p value assinged'); return() }
        if ( length(pos) != 0 ) # if there is enough data to check
        {
          if (print_pvals)
          {
            print(sprintf('%s vs %s p.val = %f , BH = %f', var_1, var_2, p_vals[pos], BH[pos] ))
          }
          if ( BH[pos] < p_val_cutoff ) # if p.value is significant enough, add asterix
          {
            # paramters
            astx_pos = max_pos
            start_pos = i
            end_pos = j
            if (manual_label == '')
            {
              lbl = find_pVal_asterix(BH[pos])
              if (!show_pval_as_asterix)
              {
                lbl = ifelse( BH[pos] < 0.009, sprintf('Q = %.2e',BH[pos]), sprintf('Q = %.4s',BH[pos]) )
                # if ( p_vals[pos] == BH[pos] )
                if ( l == 2 )  
                  lbl = gsub('Q','P',lbl)
              }
            }  else 
              lbl = manual_label
            
            # building neccecary data frame
            x_pos = c(start_pos, start_pos,end_pos,end_pos)
            y_pos = c(astx_pos, astx_pos+scale*1, astx_pos+scale*1, astx_pos)
            df2 <- data.frame(x_pos = x_pos, y_pos = y_pos )
            
            # adding the asterix to plot
            p = p + geom_path(data = df2, aes(x = x_pos, y = y_pos),size = line_width) + 
              annotate("text", x = mean( c(start_pos, end_pos) ), y = astx_pos+scale*asterix_scale_var, label = lbl, size = label_size) 
            
            # updating maximal position on plot for next time
            max_pos = max_pos + scale*2.5
            if (log_scale_flag)
              max_pos = max_pos + scale*max_pos
          }
        }
      }
    }
  }
  # creating a dataframe with the p value the BH corrected values per varibles 
  pvals_df = data.frame( var_1 = vars_1, var_2 = vars_2, p_val = p_vals, BH = BH )
  return( list(p, pvals_df) )
}

# adds significant asterix by p. value to plot, 
# by Benjamini–Hochberg procedure (not Bonferroni correncted here) and wilcox test 
# p = the plot to add on. has to be built as ggplot(data = ?), with aes in geom_boxplot etc.
# var = the variable to divide by. needs to be factor
# val = the values
# p_val_cutoff = cuoff of p.value to show differnce form. default is 0.05
# line_width = width of the line. astetic thing
# v2 - added new test type - self made permutation (ranked mean difference, 1000+ permutaion, x>= my value ,x+1/x+1) test, as Amnon recommended.
add_significant_asterix_to_plot_v3 = function(p, var, val, p_val_cutoff = 0.05, line_width = 1, print_pvals = F, log_scale_flag=F, test_type = 'perm',label_size = 8, scale_var = 6, asterix_scale_var = 1.3, manual_label = '', show_pval_as_asterix = T, fdr_flag= T)
{
  all_vars = levels(var)
  
  # calculate the number of comparisons
  l = length( all_vars )
  
  scale = sd(val, na.rm = T)/scale_var
  # if ( log_scale_flag )
  #   scale = 2^scale
  
  max_pos = max(val, na.rm = T) + scale
  # going over all comparisons between differnt variable (no repeats)
  # this repeat is to calculate p vals
  
  p_vals = c()
  vars_1 = c()
  vars_2 = c()
  is = c()
  js = c()
  for ( i in 1:l )
  {
    if ( i-1 != 0 ) # because 1:0 is 1 0 and not nothing in R
    {
      for ( j in 1:(i-1) )
      {
        var_1 = all_vars[i]
        var_2 = all_vars[j]
        
        pos1 = var == var_1 & !is.na(var)
        pos2 = var == var_2 & !is.na(var)
        
        # if there is enough data to check
        res = list()
        if ( sum(pos1)>1 & sum(pos2)>1 )
        {
          # calculate p.vlaue
          if ( test_type == 't_test') {
            ##                            big problem!!!!!!!!!!!
            # tryCatch( res = t.test(x = val[pos1], y = val[pos2], na.rm = T),
            #           error = function(e){print(sprintf('t_test error %s',e$message) )}, 
            #           finally = {res$p.value = 1})
            res = t.test(x = val[pos1], y = val[pos2], na.rm = T)
          } else if ( test_type == 'wilcox' ) {
            res = wilcox.test(x = val[pos1], y = val[pos2], na.rm = T)
          } else if ( test_type == 'perm' )
          {
            res = perm_test(x = val[pos1], y = val[pos2])
          } else {print ('test type does not fit to "t_test", "wilcox" or "perm". ')}
        } else {res$p.value = NA}
        
        p_val = res$p.value 
        if ( !is.na(p_val) ) # if a legit p value (enough data to make a comparison, etc.), update it
        {
          p_vals = c(p_vals, p_val)
          vars_1 = c(vars_1, var_1)
          vars_2 = c(vars_2, var_2)
          is = c(is, i)
          js = c(js, j)
        }
      }
    }
  }
  # calculating Benjamini–Hochberg
  BH = p.adjust(p_vals,method="BH") 
  if (! fdr_flag)
    BH = p_vals
  # going over all comparisons between differnt variable (no repeats)
  for ( i in 1:l )
  {
    if ( i-1 != 0 ) # because 1:0 is 1 0 and not nothing in R
    {
      for ( j in 1:(i-1) )
      {
        var_1 = all_vars[i]
        var_2 = all_vars[j]
        
        pos = which( is == i & js == j )
        
        if ( length(pos) >1 ) { print('ERROR! check the p value assinged'); return() }
        if ( length(pos) != 0 ) # if there is enough data to check
        {
          if (print_pvals & fdr_flag)
          {
            print(sprintf('%s vs %s p.val = %f , BH = %f', var_1, var_2, p_vals[pos], BH[pos] ))
          } else if (print_pvals)
            print(sprintf('%s vs %s p.val = %f', var_1, var_2, p_vals[pos]))
          
          if ( BH[pos] < p_val_cutoff ) # if p.value is significant enough, add asterix
          {
            # paramters
            astx_pos = max_pos
            start_pos = i
            end_pos = j
            if (manual_label == '')
            {
              lbl = find_pVal_asterix(BH[pos])
              if (!show_pval_as_asterix)
              {
                lbl = ifelse( BH[pos] < 0.009, sprintf('Q = %.2e',BH[pos]), sprintf('Q = %.4s',BH[pos]) )
                # if ( p_vals[pos] == BH[pos] )
                if ( l == 2 )  
                  lbl = gsub('Q','P',lbl)
              }
            }  else 
              lbl = manual_label
            
            # building neccecary data frame
            x_pos = c(start_pos, start_pos,end_pos,end_pos)
            y_pos = c(astx_pos, astx_pos+scale*1, astx_pos+scale*1, astx_pos)
            df2 <- data.frame(x_pos = x_pos, y_pos = y_pos )
            
            # adding the asterix to plot
            p = p + geom_path(data = df2, aes(x = x_pos, y = y_pos),size = line_width) + 
              annotate("text", x = mean( c(start_pos, end_pos) ), y = astx_pos+scale*asterix_scale_var, label = lbl, size = label_size) 
            
            # updating maximal position on plot for next time
            max_pos = max_pos + scale*2.5
            if (log_scale_flag)
              max_pos = max_pos + scale*max_pos
          }
        }
      }
    }
  }
  # creating a dataframe with the p value the BH corrected values per varibles 
  pvals_df = data.frame( var_1 = vars_1, var_2 = vars_2, p_val = p_vals, BH = BH )
  return( list(p, pvals_df) )
}

# running a permutation test, as Amnon reccomended.
# ranked mean difference, 1000+ permutaion, x>= my value ,x+1/x+1
perm_test = function(x, y, perm_num = 999)
{
  x=x[!is.na(x)]
  y=y[!is.na(y)]
  
  xy_rank = rank(c(x,y))
  x = xy_rank[1:length(x)]
  y = xy_rank[ (length(x)+1):length(xy_rank) ]
  diff = abs(mean(x)-mean(y))
  
  count=0
  for ( i in 1:perm_num  )
  {
    perm_xy = sample( c( x,y) )
    perm_x = perm_xy[1:length(x)]
    perm_y = perm_xy[ (length(x)+1):length(perm_xy) ]
    perm_diff = abs(mean(perm_x)-mean(perm_y))
    if ( perm_diff >= diff)
      count = count+1
  }
  res = list()
  res$p.value = (count+1)/(perm_num+1)
  return(res)
}

# adds significant asterix by p. value to plot, by Bonferroni correncted t/wilcox test. checks only compared to one wanted group
# p = the plot to add on. has to be built as ggplot(data = ?), with aes in geom_boxplot etc.
# var = the variable to divide by. needs to be factor
# val = the values
# group2check = the name of the group in var to check against
# p_val_cutoff = cuoff of p.value to show differnce form. default is 0.05
# line_width = width of the line. astetic thing
add_significant_asterix_to_plot_by1group = function(p, var, val, group2check, p_val_cutoff = 0.05, line_width = 1, print_pvals = F, log_scale_flag=F, test_type = 't_test')
{
  all_vars = levels(var)
  
  # calculate the number of comparisons for multiple comparisons Bonferroni correction
  l = length( all_vars )
  fdr_fix = l-1
  
  scale = sd(val, na.rm = T)/6
  # if ( log_scale_flag )
  #   scale = 2^scale
  
  max_p = -1
  max_pos = max(val, na.rm = T) + scale
  # going over all comparisons between differnt variable (no repeats)
  for ( i in 1:l )
  {
    if ( all_vars[i] != group2check )
    {
      var_1 = group2check
      var_2 = all_vars[i]
      
      pos1 = var == var_1 & !is.na(var)
      pos2 = var == var_2 & !is.na(var)
      
      # if there is enough data to check
      res = list()
      if ( sum(pos1)>1 & sum(pos2)>1 )
      {
        # calculate p.vlaue
        if ( test_type == 't_test') {
          res = t.test(x = val[pos1], y = val[pos2], na.rm = T)
        } else if ( test_type == 'wilcox' ) {
          res = wilcox.test(x = val[pos1], y = val[pos2], na.rm = T)
        } else {print ('test type does not fit to "t_test" or "wilcox')}
        if (print_pvals)
        {
          print(sprintf('%s vs %s p.val = %e', var_1, var_2, res$p.value * fdr_fix ))
          if ( res$p.value* fdr_fix < p_val_cutoff )
            if ( res$p.value* fdr_fix > max_p )
              max_p = res$p.value* fdr_fix
        }
      } else {res$p.value = 1}
      # res = wilcox.test(x = val[var == var_1], y = val[var == var_2])
      p_val = res$p.value * fdr_fix
      if ( p_val < p_val_cutoff ) # if p.value is significant enough, add asterix
      {
        # paramters
        astx_pos = max_pos
        start_pos = i
        end_pos = which(all_vars == group2check)
        lbl = find_pVal_asterix(p_val)
        
        # building neccecary data frame
        x_pos = c(start_pos, start_pos,end_pos,end_pos)
        y_pos = c(astx_pos, astx_pos+scale*1, astx_pos+scale*1, astx_pos)
        df2 <- data.frame(x_pos = x_pos, y_pos = y_pos )
        
        # adding the asterix to plot
        p = p + geom_path(data = df2, aes(x = x_pos, y = y_pos),size = line_width) + 
          annotate("text", x = mean( c(start_pos, end_pos) ), y = astx_pos+scale*1.3, label = lbl, size = 8) 
        
        # updating maximal position on plot for next time
        max_pos = max_pos + scale*2.5
        if (log_scale_flag)
          max_pos = max_pos + scale*max_pos
      }
      
    }
  }
  if (print_pvals)
    print(sprintf('maximal significant p.val = %e', max_p ))
  return(p)
}

# filter table by wanted level. written for picruat function table but shoukd work for others too.
filter_biom_table_by_level = function(taxa, level, level_sep = '__')
{
  cols = names(taxa)
  
  # get positions of metadat and wanted level
  pos = which(cols == 'Description')
  met_cols_pos = 1:pos
  
  temp_lvls = str_count(cols, level_sep)
  level_cols_pos = which(temp_lvls == level-1)
  
  pos = c(met_cols_pos, level_cols_pos)
  pos = unique(pos)
  
  return(taxa[,pos])
}

# plots boxplots of proteobacteria, firmicutes and bacteroidetes RA by check_val. check significant if wanted
pfb_RA_boxplot = function(taxa, check_val, test_type = 't_test', test_flag = TRUE)
{
  set_RA_figure = function(fig, ttl = NA)
  {
    # names = c('Israeli_Hospitalized','Israeli_Healthy','Israeli_Control','Malawi_Healthy','Venezuela_Healthy','America_Healthy','America_healthy_students')
    # new_names = c('Israeli hospitalized','Israeli healthy 1','Israeli healthy 2','Malawians','Amerindians','US' ,'US students')
    # names = c('Israeli_Hospitalized','Israeli_Healthy','Israeli_Control')
    # new_names = c('Hospitalized','Healthy 1','Healthy 2')
    # 
    fig = format_fig(fig,ttl,x_angle = 45) + theme(legend.position = 'none')
    # fig = fig + # scale_x_discrete(breaks=names, labels=new_names) + 
    #   scale_fill_manual(values = c("#3288bd","#fee08b","#f46d43"))
    return(fig)
  }
  
  test_type = 't_test'
  
  res = set_taxa_for_area_plot(taxa, file_flag = F, x_parm = check_val, taxa_cutoff = 0.01)
  plot_data = res[[1]]
  b_data = plot_data[plot_data$taxa == 'k__Bacteria.p__Bacteroidetes',]
  b_plot = ggplot() + geom_boxplot(data = b_data, aes(y = value, x = b_data[[check_val]], fill = b_data[[check_val]])) + ylab('Bacteroidetes RA') + xlab(check_val)
  b_plot = add_significant_asterix_to_plot(b_plot, b_data[[check_val]], b_data$value, print_pvals = T, test_type = test_type)
  b_plot = set_RA_figure(b_plot,'B')
  
  f_data = plot_data[plot_data$taxa == 'k__Bacteria.p__Firmicutes',]
  f_plot = ggplot() + geom_boxplot(data = f_data, aes(y = value, x = b_data[[check_val]], fill = b_data[[check_val]])) + ylab('Firmicutes RA') + xlab(check_val)
  f_plot = add_significant_asterix_to_plot(f_plot, b_data[[check_val]], f_data$value, print_pvals = T, test_type = test_type)
  f_plot = set_RA_figure(f_plot,'C')
  
  p_data = plot_data[plot_data$taxa == 'k__Bacteria.p__Proteobacteria',]
  p_data = droplevels(p_data)
  p_plot = ggplot() + geom_boxplot(data = p_data, aes(y = value, x = b_data[[check_val]], fill = b_data[[check_val]])) + ylab('Proteobacteria RA') + xlab(check_val)
  p_plot = add_significant_asterix_to_plot(p_plot, b_data[[check_val]], p_data$value, print_pvals = T, test_type = test_type)
  # p_plot = add_significant_asterix_to_plot_by1group(p_plot, p_data$Type, p_data$value, group2check = 'Israeli_Hospitalized', print_pvals = T, test_type = test_type)
  p_plot = set_RA_figure(p_plot,'A')
  
  return(list(p_plot, b_plot, f_plot))
  # return(lay_out(list(p_plot , 1, 1),
  #                list(b_plot, 2, 1),
  #                list(f_plot, 3,1)))
  
  # lay_out(list(p_plot , 1, 1),
  #         list(b_plot, 1, 2),
  #         list(f_plot, 1,3))
}
add_significant_asterix_to_plot_bothType = function(plot,var, val, p_val_cutoff = 0.05, line_width = 1, print_pvals = F, log_scale_flag=F, test_type = 't_test', label_size = 8, scale_var = 6, asterix_scale_var = 1.3, fdr = 'BH')
{
  if ( fdr == 'bonferroni' )
  {
    res = add_significant_asterix_to_plot_v2(plot, var, val, p_val_cutoff = p_val_cutoff, line_width = line_width, print_pvals = print_pvals, log_scale_flag=log_scale_flag, test_type = test_type, label_size = label_size, scale_var = scale_var, asterix_scale_var = asterix_scale_var)
  } else if ( fdr == 'BH' )
  {
    res = add_significant_asterix_to_plot_BH(plot, var, val, p_val_cutoff = p_val_cutoff, line_width = line_width, print_pvals = print_pvals, log_scale_flag=log_scale_flag, test_type = test_type, label_size = label_size, scale_var = scale_var, asterix_scale_var = asterix_scale_var)
  }
  return(res)
}

# with actinobacteria
pfba_RA_boxplot = function(taxa, check_val, test_type = 't_test', test_flag = TRUE, fdr = 'bonferroni', rm_na = T)
{
  set_RA_figure = function(fig, ttl = NA)
  {
    # names = c('Israeli_Hospitalized','Israeli_Healthy','Israeli_Control','Malawi_Healthy','Venezuela_Healthy','America_Healthy','America_healthy_students')
    # new_names = c('Israeli hospitalized','Israeli healthy 1','Israeli healthy 2','Malawians','Amerindians','US' ,'US students')
    # names = c('Israeli_Hospitalized','Israeli_Healthy','Israeli_Control')
    # new_names = c('Hospitalized','Healthy 1','Healthy 2')
    # 
    fig = format_fig(fig,ttl,x_angle = 45) + theme(legend.position = 'none')
    # fig = fig + # scale_x_discrete(breaks=names, labels=new_names) + 
    #   scale_fill_manual(values = c("#3288bd","#fee08b","#f46d43"))
    return(fig)
  }
  
  # test_type = 't_test'
  
  res = set_taxa_for_area_plot(taxa, file_flag = F, x_parm = check_val, taxa_cutoff = 0.001, rm_na = rm_na)
  plot_data = res[[1]]
  
  b_data = plot_data[plot_data$taxa == 'k__Bacteria.p__Bacteroidetes',]
  b_plot = ggplot() + geom_boxplot(data = b_data, aes(y = value, x = b_data[[check_val]], fill = b_data[[check_val]])) + ylab('Bacteroidetes RA') + xlab(check_val)
  # b_plot = add_significant_asterix_to_plot(b_plot, b_data[[check_val]], b_data$value, print_pvals = T, test_type = test_type)
  b_plot = add_significant_asterix_to_plot_bothType(b_plot, b_data[[check_val]], b_data$value, print_pvals = T, test_type = test_type, fdr =fdr)
  b_plot = b_plot[[1]]
  b_plot = set_RA_figure(b_plot,'B')
  
  f_data = plot_data[plot_data$taxa == 'k__Bacteria.p__Firmicutes',]
  f_plot = ggplot() + geom_boxplot(data = f_data, aes(y = value, x = b_data[[check_val]], fill = b_data[[check_val]])) + ylab('Firmicutes RA') + xlab(check_val)
  # f_plot = add_significant_asterix_to_plot(f_plot, b_data[[check_val]], f_data$value, print_pvals = T, test_type = test_type)
  f_plot = add_significant_asterix_to_plot_bothType(f_plot, f_data[[check_val]], f_data$value, print_pvals = T, test_type = test_type, fdr =fdr)
  f_plot = f_plot[[1]]
  f_plot = set_RA_figure(f_plot,'C')
  
  p_data = plot_data[plot_data$taxa == 'k__Bacteria.p__Proteobacteria',]
  p_data = droplevels(p_data)
  p_plot = ggplot() + geom_boxplot(data = p_data, aes(y = value, x = b_data[[check_val]], fill = b_data[[check_val]])) + ylab('Proteobacteria RA') + xlab(check_val)
  # p_plot = add_significant_asterix_to_plot(p_plot, b_data[[check_val]], p_data$value, print_pvals = T, test_type = test_type)
  p_plot = add_significant_asterix_to_plot_bothType(p_plot, p_data[[check_val]], p_data$value, print_pvals = T, test_type = test_type, fdr =fdr)
  p_plot = p_plot[[1]]
  p_plot = set_RA_figure(p_plot,'A')
  
  # browser()
  a_data = plot_data[plot_data$taxa == 'k__Bacteria.p__Actinobacteria',]
  a_plot = ggplot() + geom_boxplot(data = a_data, aes(y = value, x = a_data[[check_val]], fill = a_data[[check_val]])) + ylab('Actinobacteria RA') + xlab(check_val)
  # a_plot = add_significant_asterix_to_plot(a_plot, a_data[[check_val]], a_data$value, print_pvals = T, test_type = test_type)
  a_plot = add_significant_asterix_to_plot_bothType(a_plot, a_data[[check_val]], a_data$value, print_pvals = T, test_type = test_type, fdr =fdr)
  a_plot = a_plot[[1]]
  a_plot = set_RA_figure(a_plot,'D')
  
  return(list(p_plot, b_plot, f_plot, a_plot))
  
  # lay_out(list(p_plot , 1, 1),
  #         list(b_plot, 1, 2),
  #         list(f_plot, 1,3))
}

# plot 4 main phylum, scatter or boxplot by check_val
pfba_RA_plot = function(taxa, check_val, test_type = 'wilcox', test_flag = TRUE, fdr = 'BH', rm_na = T)
{
  set_RA_figure = function(fig, ttl = NA)
  {
    fig = format_fig(fig,ttl,x_angle = 45) + theme(legend.position = 'none')
    return(fig)
  }
  
  y_val = 'k__Bacteria.p__Bacteroidetes'; y_lb = 'Bacteroidetes'
  b_plot = plot_variable_and_stats(taxa, check_val, y_val =y_val, y_lb = y_lb, sig_ast_flag = T, print_pvals = T)
  b_plot = set_RA_figure(b_plot,'B')
  
  y_val = 'k__Bacteria.p__Firmicutes'; y_lb = 'Firmicutes'
  f_plot = plot_variable_and_stats(taxa, check_val, y_val =y_val, y_lb = y_lb, sig_ast_flag = T, print_pvals = T)
  f_plot = set_RA_figure(f_plot,'C')
  
  y_val = 'k__Bacteria.p__Proteobacteria'; y_lb = 'Proteobacteria'
  p_plot = plot_variable_and_stats(taxa, check_val, y_val =y_val, y_lb = y_lb, sig_ast_flag = T, print_pvals = T)
  p_plot = set_RA_figure(p_plot,'A')
  
  y_val = 'k__Bacteria.p__Actinobacteria'; y_lb = 'Actinobacteria'
  a_plot = plot_variable_and_stats(taxa, check_val, y_val =y_val, y_lb = y_lb, sig_ast_flag = T, print_pvals = T)
  a_plot = set_RA_figure(a_plot,'D')
  
  return(list(p_plot, b_plot, f_plot, a_plot))
  
  # lay_out(list(p_plot , 1, 1),
  #         list(b_plot, 1, 2),
  #         list(f_plot, 1,3))
}
# with actinobacteria
fb_ratio_p_RA_boxplot = function(taxa, check_val, test_type = 't_test', test_flag = TRUE, fdr = 'bonferroni')
{
  set_RA_figure = function(fig, ttl = NA)
  {
    # names = c('Israeli_Hospitalized','Israeli_Healthy','Israeli_Control','Malawi_Healthy','Venezuela_Healthy','America_Healthy','America_healthy_students')
    # new_names = c('Israeli hospitalized','Israeli healthy 1','Israeli healthy 2','Malawians','Amerindians','US' ,'US students')
    # names = c('Israeli_Hospitalized','Israeli_Healthy','Israeli_Control')
    # new_names = c('Hospitalized','Healthy 1','Healthy 2')
    # 
    fig = format_fig(fig,ttl,x_angle = 45) + theme(legend.position = 'none')
    # fig = fig + # scale_x_discrete(breaks=names, labels=new_names) + 
    #   scale_fill_manual(values = c("#3288bd","#fee08b","#f46d43"))
    return(fig)
  }
  
  # test_type = 't_test'
  
  res = set_taxa_for_area_plot(taxa, file_flag = F, x_parm = check_val, taxa_cutoff = 0.001)
  plot_data = res[[1]]
  
  b_data = plot_data[plot_data$taxa == 'k__Bacteria.p__Bacteroidetes',]
  b_plot = ggplot() + geom_boxplot(data = b_data, aes(y = value, x = b_data[[check_val]], fill = b_data[[check_val]])) + ylab('Bacteroidetes RA') + xlab(check_val)
  # b_plot = add_significant_asterix_to_plot(b_plot, b_data[[check_val]], b_data$value, print_pvals = T, test_type = test_type)
  b_plot = add_significant_asterix_to_plot_bothType(b_plot, b_data[[check_val]], b_data$value, print_pvals = T, test_type = test_type, fdr =fdr)
  b_plot = b_plot[[1]]
  b_plot = set_RA_figure(b_plot,'C')
  
  f_data = plot_data[plot_data$taxa == 'k__Bacteria.p__Firmicutes',]
  f_plot = ggplot() + geom_boxplot(data = f_data, aes(y = value, x = b_data[[check_val]], fill = b_data[[check_val]])) + ylab('Firmicutes RA') + xlab(check_val)
  # f_plot = add_significant_asterix_to_plot(f_plot, b_data[[check_val]], f_data$value, print_pvals = T, test_type = test_type)
  f_plot = add_significant_asterix_to_plot_bothType(f_plot, f_data[[check_val]], f_data$value, print_pvals = T, test_type = test_type, fdr =fdr)
  f_plot = f_plot[[1]]
  f_plot = set_RA_figure(f_plot,'B')
  
  fb_data = f_data
  fb_data$value = f_data$value / b_data$value
  fb_plot = ggplot() + geom_boxplot(data = fb_data, aes(y = value, x = fb_data[[check_val]], fill = fb_data[[check_val]])) + ylab('Firmicutes/Bacteroidetes RA') + xlab(check_val)
  # f_plot = add_significant_asterix_to_plot(f_plot, b_data[[check_val]], f_data$value, print_pvals = T, test_type = test_type)
  fb_plot = add_significant_asterix_to_plot_bothType(fb_plot, fb_data[[check_val]], fb_data$value, print_pvals = T, test_type = test_type, fdr =fdr)
  fb_plot = fb_plot[[1]]
  fb_plot = set_RA_figure(fb_plot,'D')
  
  p_data = plot_data[plot_data$taxa == 'k__Bacteria.p__Proteobacteria',]
  p_data = droplevels(p_data)
  p_plot = ggplot() + geom_boxplot(data = p_data, aes(y = value, x = b_data[[check_val]], fill = b_data[[check_val]])) + ylab('Proteobacteria RA') + xlab(check_val)
  # p_plot = add_significant_asterix_to_plot(p_plot, b_data[[check_val]], p_data$value, print_pvals = T, test_type = test_type)
  p_plot = add_significant_asterix_to_plot_bothType(p_plot, p_data[[check_val]], p_data$value, print_pvals = T, test_type = test_type, fdr =fdr)
  p_plot = p_plot[[1]]
  p_plot = set_RA_figure(p_plot,'A')
  
  return(list(p_plot, b_plot, f_plot, fb_plot, p_plot))
  
  # lay_out(list(p_plot , 1, 1),
  #         list(b_plot, 1, 2),
  #         list(f_plot, 1,3))
}
# return level of taxa id
check_id_level = function( id )
{
  return(str_count(id,'__'))
}

rainbow_colour_gradient = function(ttl = NA)
{
  # unqiue colours
  library(RColorBrewer)
  n <- 40
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  # pie(rep(1,n), col=sample(col_vector, n))
  col_vector = col_vector[1:21]
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  if (is.na(ttl)) {
    sc <- scale_colour_gradientn(colours = myPalette(100))
  } else {
    sc <- scale_colour_gradientn(colours = myPalette(100), name = ttl)
  }
  return(sc)
}

# creates a ggplot scale_colour_discrete() object with rainbow colors and wanted title/name
rainbow_colour_discrete = function(ttl = NA)
{
  # # unqiue colours
  # library(RColorBrewer)
  # n <- 40
  # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  # # pie(rep(1,n), col=sample(col_vector, n))
  # col_vector = col_vector[1:21]
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  if (is.na(ttl)) {
    sc <- scale_colour_discrete(colours = myPalette(100))
  } else {
    sc <- scale_colour_discrete(colours = myPalette(100), name = ttl)
  }
  return(sc)
}

# creates a ggplot scale_colour_discrete() object with rainbow colors and wanted title/name
unique_color_manual = function(n, ttl = NA)
{
  # unqiue colours
  library(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  # pie(rep(1,n), col=sample(col_vector, n))
  col_vector = col_vector[1:n]
  if (is.na(ttl)) {
    sc <- scale_color_manual(values=col_vector)
  } else {
    sc <- scale_color_manual(values=col_vector, name = ttl)
  }
  return( sc ) 
}


# plots the highest abundance taxa, differentiated by wanted parameter.
# gets taxa dataframe, number of top taxa to include, and the wanted differentiation parameter.
# returns the plot, the taxa names and the used new dataframe with the data.
highest_aboundance_taxa_plot = function(taxa, number, col_prm, median_flag = F, picrust_flag = F)
{
  if ( picrust_flag )
  {
    taxa_pos = which( grepl(pattern = '__.*__',x = names(taxa)) )
  } else
  {
    taxa_pos = which( grepl(pattern = 'g__',x = names(taxa)) & ! grepl(pattern = 's__',x = names(taxa)) )
  }
  sums = colSums(taxa[,taxa_pos])
  if (median_flag)
    sums = apply(taxa[,taxa_pos], 2, median)
  highest_15_pos = rev(order(sums))[1:number]
  highest_names = names(taxa)[ taxa_pos[highest_15_pos] ]
  highest_df = taxa[highest_names]
  
  highest_df[[col_prm]] = taxa[[col_prm]]
  highest_df = melt(data = highest_df, id.vars = col_prm)
  
  highest_df$variable = factor(highest_df$variable, levels = rev(levels(highest_df$variable)) )
  highest_df = highest_df[ !is.na(highest_df[[col_prm]] ), ]
  plot = ggplot(highest_df, aes(x = variable, y = value)) +  geom_boxplot(aes(fill = highest_df[[col_prm]]), alpha = 0.6) + 
    coord_flip()  + geom_point(position=position_jitterdodge(), aes(colour = highest_df[[col_prm]])) + 
    ylab('Abundance') + xlab('') # + scale_y_log10()
  if ( median_flag )
    plot = plot + annotate("text", x=0.1, y=Inf, hjust=0, vjust=1,label = 'median') 
  return(list(plot, highest_names, highest_df))
}


# gets a data frame with beta diversity between samples, all from same patient. 
# finds the time point the is on average the most different from the others.
calc_furthest_time = function(bDiv_df)
{
  lvls = levels(bDiv_df$time1)
  
  max_dist = 0
  far_lvl = NA
  for ( lvl in lvls )
  {
    vals = bDiv_df$beta_dsitance[ bDiv_df$time1 == lvl | bDiv_df$time2 == lvl ]
    lvl_dist = mean( vals, na.rm = T )
    if ( lvl_dist > max_dist )
    {
      max_dist = lvl_dist
      far_lvl = lvl
    }
  }
  pos = which( far_lvl == lvls)
  return(as.numeric(pos))
}

# make a heatmap of the beta diversity between all samples, by variable (built for panient samples over time)
bDiv_samples_heatmap = function(b_file, taxa, check_val = 'Stool_number', weight = 'unweighted', na_str = c('no_data','_','NA','unknown', 'other'), Dx_flag = T )
{
  # read beta div file
  bDiv = read.table(b_file,sep="\t", header=T, na.strings = na_str, row.names = 'X')
  
  lvls = levels( as.factor(taxa[[check_val]]) )
  lvls = sort(lvls)
  if ( is.numeric(taxa[[check_val]]) )
  {
    lvls = unique(taxa[[check_val]])
    lvls = sort(lvls)
    lvls = as.character( lvls )
  }
  
  vals=  c()
  time1 = c()
  time2 = c()
  for ( i in 1:length(lvls) )
  {
    id1 = as.character( taxa$SampleID[ taxa[[check_val]] == lvls[[i]] ] )
    for ( j in 1:length(lvls)  )
    {
      
      if ( i < j )
      {
        
        id2 = as.character( taxa$SampleID[ taxa[[check_val]] == lvls[[j]] ] )
        
        val = bDiv[id1, id2]
        
        vals = c(vals, val)
      } 
      # adding marker to Dx change in the i==j diagonal
      if ( i==j ) 
      {
        if ( Dx_flag ) 
        {
          if ( temp_taxa$Dx_PGA[ which(taxa[[check_val]] == lvls[[j]]) ] == 'CD_active' ) 
          {
            vals = c(vals, 0)
          } else {
            vals = c(vals, NA)
          }
        } else 
          vals = c(vals, NA)
      } 
      if ( i>j ) 
      {
        vals = c(vals, NA)
        
        # temp = temp_taxa$dysbiosis_index[ taxa[[check_val]] == lvls[[i]] ] - temp_taxa$dysbiosis_index[ taxa[[check_val]] == lvls[[j]] ]
        # vals = c(vals, temp  )
        
      } 
      time1 = c(time1, lvls[[i]])
      time2 = c(time2, lvls[[j]])
      
    }
  }
  bDiv_df = data.frame(time1 = time1, time2 = time2, beta_dsitance = vals)
  bDiv_df$patient = taxa$patient_No[1]
  pos = calc_furthest_time(bDiv_df)
  
  plot = ggplot(bDiv_df, aes( x = time1, y=time2, fill = beta_dsitance) ) + geom_tile() + 
    scale_y_discrete(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0)) +
    scale_fill_gradientn(colours = c("red","white","blue"), 
                         guide = "colorbar", 
                         name = 'beta diversity\ndistance\n' , limits = c(0,1)) + 
    annotate("text", x=pos, y=pos, label = 'far') 
  
  
  return(list( plot, bDiv_df ))
}

# calculate the average beta diversity between adjusent sampels 
calc_close_bDiv = function( bDiv_df )
{
  bDiv_df = bDiv_df[ !is.na(bDiv_df$beta_dsitance), ]
  pos = which( abs(as.numeric(bDiv_df$time1)-as.numeric(bDiv_df$time2) ) == 1 )
  return( mean(bDiv_df$beta_dsitance[pos]) )
}

# calculate the average beta diversity between adjusent sampels 
calc_close_bDiv_med = function( bDiv_df )
{
  bDiv_df = bDiv_df[ !is.na(bDiv_df$beta_dsitance), ]
  pos = which( abs(as.numeric(bDiv_df$time1)-as.numeric(bDiv_df$time2) ) == 1 )
  return( median(bDiv_df$beta_dsitance[pos]) )
}

# calculate the average beta diversity between adjusent sampels 
calc_close_bDiv_min = function( bDiv_df )
{
  bDiv_df = bDiv_df[ !is.na(bDiv_df$beta_dsitance), ]
  pos = which( abs(as.numeric(bDiv_df$time1)-as.numeric(bDiv_df$time2) ) == 1 )
  return( min(bDiv_df$beta_dsitance[pos]) )
}

# calculate the average beta diversity between sampels 
calc_all_bDiv = function( bDiv_df )
{
  bDiv_df = bDiv_df[ !is.na(bDiv_df$beta_dsitance), ]
  pos = which( as.numeric(bDiv_df$time1) < as.numeric(bDiv_df$time2 ) )
  # return( min(bDiv_df$beta_dsitance[pos]) )
  return( mean(bDiv_df$beta_dsitance[pos]) )
}

# calculate the average beta diversity between sampels 
calc_first_last_bDiv = function( bDiv_df )
{
  bDiv_df = bDiv_df[ !is.na(bDiv_df$beta_dsitance), ]
  pos = which( as.numeric(bDiv_df$time1) == 1 & (as.numeric(bDiv_df$time2) == max(as.numeric(bDiv_df$time2)) )  )
  return( bDiv_df$beta_dsitance[pos] )
}

take_last_samples = function(temp_taxa, samples_num = 3)
{
  temp_taxa = temp_taxa[temp_taxa$Stool_number > max( temp_taxa$Stool_number ) -samples_num, ]
  return(temp_taxa)
}

take_first_samples = function(temp_taxa, samples_num = 3)
{
  temp_taxa = temp_taxa[temp_taxa$Stool_number <= samples_num, ]
  return(temp_taxa)
}

piecrust_RA_boxplot = function(taxa, check_val, test_type = 't_tets', test_flag = TRUE)
{
  phyls = c('Cellular_Processes','Environmental_Information_Processing',
            'Genetic_Information_Processing',
            'Human_Diseases','Metabolism',
            'Organismal_Systems','Unclassified')
  
  RA_res = list()
  temp_taxa = taxa[ !is.na(taxa[[check_val]]), ]
  for ( i in 1:length(phyls) )
  {
    var = temp_taxa[[check_val]]
    val = temp_taxa[[ phyls[i] ]]
    
    p = ggplot(temp_taxa) + 
      geom_boxplot(aes_string(x=check_val, y=phyls[i] )) + 
      theme_bw() + xlab(check_val) +  ylab(phyls[i])
    
    if ( test_flag )
    {
      res = add_significant_asterix_to_plot_v2(p = p, var = var, val = val,test_type = test_type)
      p=res[[1]]
    }
    RA_res[[i]] = p
  }
  return(RA_res)
}

test_2_vals = function(val, pos1, pos2, test_type = 't_test')
{
  if ( test_type == 't_test') {
    res = t.test(x = val[pos1], y = val[pos2], na.rm = T)
  } else if ( test_type == 'wilcox' ) {
    res = wilcox.test(x = val[pos1], y = val[pos2], na.rm = T)
  } else if (test_type == 'perm' ) {
    res = perm_test(x=val[pos1], y = val[pos2])
  } else {print ('test type does not fit to "t_test", "wilcox" or "perm". ')}
  return(res)
}

add_asterix_after_calculations = function(p, p_val, start_pos, end_pos, astx_pos, lbl, scale, line_width )
{
  lbl = find_pVal_asterix(p_val)
  
  # building neccecary data frame
  x_pos = c(start_pos, start_pos,end_pos,end_pos)
  y_pos = c(astx_pos, astx_pos+scale*1, astx_pos+scale*1, astx_pos)
  df2 <- data.frame(x_pos = x_pos, y_pos = y_pos )
  
  # adding the asterix to plot
  if  ( lbl == 'ns' )
  {
    p = p + geom_path(data = df2, aes(x = x_pos, y = y_pos),size = line_width) + 
      annotate("text", x = mean( c(start_pos, end_pos) ), y = astx_pos+scale*3, label = lbl, size = 6) 
  } else
  {
    p = p + geom_path(data = df2, aes(x = x_pos, y = y_pos),size = line_width) + 
      annotate("text", x = mean( c(start_pos, end_pos) ), y = astx_pos+scale*1.3, label = lbl, size = 8) 
  }
  return(p)
}

boxplot_2_vars_grouped_add_asterix = function( p, taxa, x_Val, y_val, fill_val, test_type = 't_test', fdr= 3)
{
  taxa[[x_Val]] = as.factor(taxa[[x_Val]])
  taxa[[fill_val]] = as.factor(taxa[[fill_val]])
  
  line_width = 1

  val = taxa[[y_val]]
  scale = sd(val, na.rm = T)/6
  
  all_vars = levels( taxa[[x_Val]] )
  all_fills = levels( taxa[[fill_val]] )

  ## within group 1
  pos1 = taxa[[x_Val]] == all_vars[1] & taxa[[fill_val]] == all_fills[1]
  pos2 = taxa[[x_Val]] == all_vars[1] & taxa[[fill_val]] == all_fills[2]
  
  res = test_2_vals(val, pos1, pos2, test_type = test_type)
  
  p_val = res$p.value * fdr
  # paramters
  astx_pos = max(val) + scale
  start_pos = 0.8
  end_pos = 1.2
  
  p = add_asterix_after_calculations(p, p_val, start_pos, end_pos, astx_pos, lbl, scale, line_width )
  
  ## within group 2
  pos1 = taxa[[x_Val]] == all_vars[2] & taxa[[fill_val]] == all_fills[1]
  pos2 = taxa[[x_Val]] == all_vars[2] & taxa[[fill_val]] == all_fills[2]
  
  # pos1 = taxa[[x_Val]] == all_vars[1]
  # pos2 = taxa[[x_Val]] == all_vars[2] 
  
  res = test_2_vals(val, pos1, pos2, test_type = test_type)
  
  p_val = res$p.value * fdr
  # paramters
  astx_pos = max(val) + scale
  start_pos = 1.8
  end_pos = 2.2
  
  p = add_asterix_after_calculations(p, p_val, start_pos, end_pos, astx_pos, lbl, scale, line_width )
  
  ## between groups
  pos1 = taxa[[x_Val]] == all_vars[1]
  pos2 = taxa[[x_Val]] == all_vars[2]
  
  res = test_2_vals(val, pos1, pos2, test_type = test_type)
  
  p_val = res$p.value * fdr
  # paramters
  astx_pos = max(val) + scale*5
  start_pos = 1
  end_pos = 2
  
  p = add_asterix_after_calculations(p, p_val, start_pos, end_pos, astx_pos, lbl, scale, line_width )
  
  return(p)
}

boxplot_2_vars_grouped = function(taxa, x_Val, y_val, fill_val, test_type = 't_test', fdr=3, add_asterix_flag=T, add_x_n_flag = F)
{
  taxa2=taxa
  taxa2 = taxa2[!is.na(taxa2[[x_Val]]), ]
  taxa2 = taxa2[!is.na(taxa2[[fill_val]]), ]
  if (add_x_n_flag)
  {
    temp =  taxa2[[x_Val]]
    for ( i in 1:length(taxa2[[x_Val]] ) )
    {
      taxa2[[x_Val]][i] = sprintf('%s (n=%d)', temp[i], sum( temp == temp[i], na.rm = T ) )
      if (is.na(temp[i]))
        taxa2[[x_Val]][i] = sprintf('%s (n=%d)', temp[i], sum( is.na(temp) ) )
    }
  }
  taxa2[[x_Val]] = as.factor(taxa2[[x_Val]])
  taxa2[[fill_val]] = as.factor(taxa2[[fill_val]])
  p = ggplot(taxa2) + geom_boxplot(aes(x=taxa2[[x_Val]], y=taxa2[[y_val]], fill = taxa2[[fill_val]])) + 
    xlab(x_Val) + ylab(y_val) + labs(fill = fill_val) + theme_bw()
  if (add_asterix_flag)
    p = boxplot_2_vars_grouped_add_asterix( p, taxa2, x_Val, y_val, fill_val, test_type = test_type,fdr=fdr)
  return(p)
}



# calculate the wanted variable inside every patient with enough samples, and return a dataframe with needed data
in_patient_var_get = function(taxa, check_var = 'bDiv', samples_number = 2, remission_samples_flag = T, last_sample_flag = T, b_file = NA, weight = 'unweighted')
{
  dists = c()
  Dxs = c()
  patients = c()
  for ( i in 1:length(levels(taxa$patient_No)) )
  {
    pn = levels(taxa$patient_No)[[i]]
    temp_taxa = taxa[taxa$patient_No == pn & !is.na(taxa$patient_No) & !taxa$full_Dx=='Control', ]
    
    if ( remission_samples_flag )
    {
      # take only samples in remission (not active)
      temp_taxa = temp_taxa[temp_taxa$Dx_PGA == 'CD_remission',]
    }
    if ( last_sample_flag )
    {
      # take last few samples
      temp_taxa = take_last_samples(temp_taxa, samples_number)
    }
    
    # take first few samples
    # temp_taxa = take_first_samples(temp_taxa,samples_number)
    
    if (dim(temp_taxa)[1] < samples_number )
    {
      next()
    }
    # browser()
    if ( check_var == 'bDiv') # if beta diversity use needed function, else just use average
    {
      res = bDiv_samples_heatmap(b_file, temp_taxa, weight = weight, Dx_flag = F)
      dist = calc_all_bDiv( res[[2]] )
    } else if ( is.numeric( temp_taxa[[check_var]] ) )
    {
      dist = mean(temp_taxa[[check_var]], na.rm = T)
    } else
    {
      dist = temp_taxa[[check_var]][1]
    }
    dists = c(dists, dist)
    
    temp = taxa$Group[taxa$patient_No == pn]
    Dxs = c(Dxs, as.character(temp[[1]])  )
    patients = c(patients, pn)
  }
  temp_df = data.frame(dists = dists, Dx = Dxs, patient = patients )
  return(temp_df)
}

# calculate the wanted variable inside every patient with enough samples, and return a dataframe with needed data
in_patient_var_get_v2 = function(taxa, check_var = 'bDiv', samples_number = 2, remission_samples_flag = T, last_sample_flag = T, b_file = NA, weight = 'unweighted', patinet_var = 'patient_No')
{
  dists = c()
  Dxs = c()
  patients = c()
  for ( i in 1:length(levels(taxa[[patinet_var]])) )
  {
    pn = levels(taxa[[patinet_var]])[[i]]
    temp_taxa = taxa[taxa[[patinet_var]] == pn & !is.na(taxa[[patinet_var]]) & !taxa$full_Dx=='Control', ]
    
    if ( remission_samples_flag )
    {
      # take only samples in remission (not active)
      temp_taxa = temp_taxa[temp_taxa$Dx_PGA == 'CD_remission',]
    }
    if ( last_sample_flag )
    {
      # take last few samples
      temp_taxa = take_last_samples(temp_taxa, samples_number)
    }
    
    # take first few samples
    # temp_taxa = take_first_samples(temp_taxa,samples_number)
    
    if (dim(temp_taxa)[1] < samples_number )
    {
      next()
    }
    # browser()
    if ( check_var == 'bDiv') # if beta diversity use needed function, else just use average
    {
      res = bDiv_samples_heatmap(b_file, temp_taxa, weight = weight, Dx_flag = F)
      dist = calc_all_bDiv( res[[2]] )
    } else if ( is.numeric( temp_taxa[[check_var]] ) )
    {
      dist = mean(temp_taxa[[check_var]], na.rm = T)
    } else
    {
      dist = temp_taxa[[check_var]][1]
    }
    dists = c(dists, dist)
    
    temp = taxa$Group[taxa[[patinet_var]] == pn]
    Dxs = c(Dxs, as.character(temp[[1]])  )
    patients = c(patients, pn)
  }
  temp_df = data.frame(dists = dists, Dx = Dxs, patient = patients )
  return(temp_df)
}
calc_flare_index = function(taxa)
{
  up_vals = taxa$k__Bacteria.p__Firmicutes.c__Bacilli.o__Gemellales.f__Gemellaceae +
    taxa$k__Bacteria.p__Tenericutes.c__Mollicutes.o__Mycoplasmatales.f__Mycoplasmataceae
  dowm_vals = taxa$k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__S24.7 +
    taxa$k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Lachnospiraceae.g__Coprococcus # + 
  # taxa$k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Christensenellaceae
  
  flare_index = up_vals / (dowm_vals + min(dowm_vals[dowm_vals!=0]))
  return(flare_index)
}

taxa_full_name_to_last_level = function(tax)
{
  name = gsub(pattern = '.*[a-z]__',replacement = '',x = as.character( tax ) )
  temp = strsplit(x = as.character(tax),split = '__')
  level = strsplit( x = temp[[1]][length(temp[[1]]) - 1 ], split = '[.]')[[1]][2]
  name = sprintf('%s__%s',level, name)
  return(name)
}

taxa_full_name_to_last_level_many = function(taxs)
{
  taxs = as.character(taxs)
  res = c()
  for (i in 1:length( taxs ) )
    if ( grepl( pattern = 'k__',x = taxs[i] ) )
    {
      res = c(res, taxa_full_name_to_last_level(taxs[i]) )
    } else
      res = c(res, taxs[i])
  return(res)
}

# 
# get minimal value in the dataset,  
# minimal value is the minimal taxonomic RA in the dataset
get_taxa_min_val = function( taxa )
{
  poses = grep(pattern = 'k__',x = names(taxa))
  
  temp =  taxa[,poses]
  temp=melt(data = temp)
  min_taxa_abd = min(temp$value[temp$value!=0])
  
  return(min_taxa_abd)
}

# change zeros in taxa aboundance to minimal value. 
# the minimal value is the minimal RA in the dataset
set_taxa_zeros_as_minVal = function( taxa )
{
  
  poses = grep(pattern = 'k__',x = names(taxa))

  min_taxa_abd = get_taxa_min_val( taxa )
  
  for (pos in poses)
  {
    taxa[,pos][ taxa[,pos]==0 ] = min_taxa_abd
  }
  return(taxa)
}


# filteres samples with NA in var, and orders both by row.name
set_dist_map_by_var = function( dist, map, var)
{
  wanted_samples_id = row.names(map)[ !is.na(map[[var]]) ]
  # set fit for how R reads names
  row.names(dist) = make.names(row.names(dist))
  dist_wanted_samples_id = make.names(wanted_samples_id)
  
  dist = dist[dist_wanted_samples_id, dist_wanted_samples_id]
  map = map[wanted_samples_id, ]
  
  # order by same order (just in case)
  dist = dist[order(row.names(dist)), order(names(dist))]
  map = map[order(row.names(map)),]
  
  return(list(dist, map))
}

set_dist_map_by_vars = function( dist, map, vars)
{
  temp = is.na(map[,vars])
  wanted_samples_id = row.names(map)[as.numeric(rowSums(temp)) == 0]
  
  # wanted_samples_id = row.names(map)[ !is.na(map[[var]]) ]
  # set fit for how R reads names
  row.names(dist) = make.names(row.names(dist))
  dist_wanted_samples_id = make.names(wanted_samples_id)
  
  dist = dist[dist_wanted_samples_id, dist_wanted_samples_id]
  map = map[wanted_samples_id, ]
  
  # order by same order (just in case)
  dist = dist[order(row.names(dist)), order(names(dist))]
  map = map[order(row.names(map)),]
  
  return(list(dist, map))
}
# for categories with only one sample group all to "other"
# if there is only one change to NA.
remove_character_list_outliers = function(v)
{
  if ( class(v) == 'numeric' | class(v) == 'integer' )
    return(v)
  if ( all( is.na(v) ) )
    return(v)
  for ( i in 1:length(v) )
  {
    if ( !is.na(v[i]) )
    {
      if ( sum( v==v[i], na.rm = T ) == 1  )
      {
        v[i] = 'other'
      }
    }
  }
  if ( sum(v=='other', na.rm = T) == 1 )
    v[v=='other'] = NA
  return(v)
}

remove_outliers_from_map = function(map,vars_2_check)
{
  vars_2_check = vars_2_check[!vars_2_check %in% c('pn_ID','pn_id') ]
  for ( i in 1:length(vars_2_check) )
  {
    map[[vars_2_check[i]]] = remove_character_list_outliers( map[[vars_2_check[i]]] )
  }
  return(map)
}

# goinf over a list of parametres in the map and checking this is more than 1 option for each of them.
# retunrs a list of "good" parameters and prints the names of the "bad" ones
check_vars_valid = function(map, vars_2_check, min_notNA_num = 0)
{
  if (length(vars_2_check) ==0)
    return(vars_2_check)
  for ( i in length(vars_2_check):1 )
  {
    temp = map[[vars_2_check[i]]]
    temp = temp[!is.na(temp)]
    l = length( unique(temp) )
    if ( l < 2 )
    {
      print( sprintf('rm %s, not valid, less than 2 unique values',vars_2_check[i]) )
      vars_2_check = vars_2_check[-i]
    } else if ( sum(!is.na(temp)) <= min_notNA_num )
    {
      print( sprintf('rm %s, not valid, %s samples not NA',vars_2_check[i], sum(!is.na(temp))) )
      vars_2_check = vars_2_check[-i]
    }
  }
  return(vars_2_check)
}

run_adonis2 = function(vars_2_check, dist, map)
{
  R2 = vector(mode = 'numeric',length = length(vars_2_check))
  pval = vector(mode = 'numeric',length = length(vars_2_check))
  for ( i in 1:length(vars_2_check) )
  {
    vr = vars_2_check[i]
    res = set_dist_map_by_var(dist, map, vr)
    dist_f = res[[1]]
    map_f = res[[2]]
    grp = map_f[[vr]]
    permanova = vegan::adonis2( dist_f~grp, permutations=999)
    # permanova = vegan::adonis2( dist_f~group, data = map_f, permutations=999)
    
    pval[i] = permanova$`Pr(>F)`[1]
    R2[i] = permanova$R2[1]
  }
  return(list(pval, R2))
}

check_var_NA = function(vars_2_check, map)
{
  var_num = vector(mode = 'character',length = length(vars_2_check))
  NA_num = vector(mode = 'numeric',length = length(vars_2_check))
  samp_num = vector(mode = 'numeric',length = length(vars_2_check))
  for ( i in 1:length(vars_2_check))
  {
    var_num[i] = as.character( sum(!is.na(( unique(map[[vars_2_check[i]]]) )) ) )
    if (is.numeric(map[[vars_2_check[i]]] ))
      var_num[i] = 'Numeric'
    NA_num[i] = sum(is.na(map[[vars_2_check[i]]]))
    samp_num[i] = sum(!is.na(map[[vars_2_check[i]]]))
  }
  return(list(var_num, NA_num,samp_num))
}

run_permanova = function(map, dist, vars_2_check)
{
  map = map[order(row.names(map)),]
  dist = dist[order(row.names(dist)), order(names(dist))]
  
  map = droplevels(map)
  
  vars_2_check = check_vars_valid(map, vars_2_check)
  
  # dist = dist[make.names(row.names(dist)) %in% make.names(row.names(map)), 
  #             make.names(row.names(dist)) %in% make.names(row.names(map))]
  # 
  source('/pita/users/tzipi/projects/16S/haya_metaanlysis/sep1/adonis_script.R')
  
  res = check_var_NA(vars_2_check, map)
  var_num = res[[1]]; NA_num = res[[2]]
  
  perm_df = data.frame(var = vars_2_check, R2 = R2, pval = pval, category_num = var_num, NA_num = NA_num)
  perm_df = perm_df[rev(order(perm_df$R2)),]
  
  return(perm_df)
}


add_manual_significant_asterix = function(g, val, lbl, start_pos = 1, end_pos = 2, scale_var = 1.3, 
                                          astx_pos = -1,asterix_scale_var = 1.3, label_size = 8, label_num_to_ast_flag = F )
{
  scale = sd(val,na.rm = T ) * scale_var
  if ( astx_pos == -1 )
    astx_pos = max(val) + scale 
  
  if ( label_num_to_ast_flag )
    if ( is.numeric(lbl))
      lbl = find_pVal_asterix(lbl)
    else
      print('label has to be numeric to be converted to an asterix.')
  
  x_pos = c(start_pos, start_pos,end_pos,end_pos)
  y_pos = c(astx_pos, astx_pos+scale*1, astx_pos+scale*1, astx_pos)
  df2 <- data.frame(x_pos = x_pos, y_pos = y_pos )
  
  # adding the asterix to plot
  g = g + geom_path(data = df2, aes(x = x_pos, y = y_pos),size = 1) + 
    annotate("text", x = mean( c(start_pos, end_pos) ), y = astx_pos+scale*asterix_scale_var, label = lbl, size = label_size) 
  return(g)
}

test_multiple_BH = function(var, val, test_type = 'wilcox')
{
  all_vars = levels(as.factor(var))
  
  # calculate the number of comparisons
  l = length( all_vars )
  
  # going over all comparisons between differnt variable (no repeats)
  # this repeat is to calculate p vals
  
  p_vals = c()
  vars_1 = c()
  vars_2 = c()
  is = c()
  js = c()
  for ( i in 1:l )
  {
    if ( i-1 != 0 ) # because 1:0 is 1 0 and not nothing in R
    {
      for ( j in 1:(i-1) )
      {
        var_1 = all_vars[i]
        var_2 = all_vars[j]
        
        pos1 = var == var_1 & !is.na(var)
        pos2 = var == var_2 & !is.na(var)
        
        # if there is enough data to check
        res = list()
        if ( sum(pos1)>1 & sum(pos2)>1 )
        {
          # calculate p.vlaue
          if ( test_type == 't_test') {
            ##                            big problem!!!!!!!!!!!
            # tryCatch( res = t.test(x = val[pos1], y = val[pos2], na.rm = T),
            #           error = function(e){print(sprintf('t_test error %s',e$message) )}, 
            #           finally = {res$p.value = 1})
            res = t.test(x = val[pos1], y = val[pos2], na.rm = T)
          } else if ( test_type == 'wilcox' ) {
            res = wilcox.test(x = val[pos1], y = val[pos2], na.rm = T)
          } else {print ('test type does not fit to "t_test" or "wilcox')}
        } else {res$p.value = NA}
        
        p_val = res$p.value 
        if ( !is.na(p_val) ) # if a legit p value (enough data to make a comparison, etc.), update it
        {
          p_vals = c(p_vals, p_val)
          vars_1 = c(vars_1, var_1)
          vars_2 = c(vars_2, var_2)
          is = c(is, i)
          js = c(js, j)
        }
      }
    }
  }
  # calculating Benjamini–Hochberg
  BH = p.adjust(p_vals,method="BH") 
  # creating a dataframe with the p value the BH corrected values per varibles 
  pvals_df = data.frame( var_1 = vars_1, var_2 = vars_2, p_val = p_vals, BH = BH )
  return(pvals_df)
}

