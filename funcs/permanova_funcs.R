
## fit map and dist dataframes (same samples, same order) and and filter wanted variable to check to remove 
## variables with to much NAs or too little unique options.
permanova_organize_data = function(map, dist, vars_2_check, min_notNA_num = 3)
{
  map = map[row.names(map) %in% row.names(dist),]
  map = droplevels(map)
  map = map[order(row.names(map)),]
  dist = dist[make.names(row.names(map)), make.names(row.names(map))]
  
  vars_2_check_here = check_vars_valid(map, vars_2_check, min_notNA_num)
  return(list(map, dist, vars_2_check_here))
}

## running permanova on each of the variables in vars_2_check_here. using map and dist data.
## controlling for extra_vars (if not c()) and stratifying by strata_var (if not NA)
permanova_run = function( map, dist, vars_2_check_here, extra_vars = c(), strata_var = NA  )
{
  R2 = vector(mode = 'numeric',length = length(vars_2_check_here))
  pval = vector(mode = 'numeric',length = length(vars_2_check_here))
  for ( i in 1:length(vars_2_check_here) )
  {
    var = vars_2_check_here[i]
    res = set_dist_map_by_vars(dist, map, c(var, extra_vars))
    dist_f = res[[1]]
    map_f = res[[2]]
    group = map_f[[var]]
    
    print(var)
    
    # create a formula accounting first for the extra vars (such as age, gender) and than for the interesting var.
    extra_vars2 = extra_vars[!extra_vars %in% var]
    prm_vars = c(extra_vars2,var)
    prm_vars = check_vars_valid(map_f, prm_vars, min_notNA_num = 3)
    frm = as.formula( sprintf('dist_f ~  %s',paste(prm_vars,collapse = ' + ') ) )
    
    if( !is.na(strata_var) )
      permanova <- vegan::adonis( frm , data = map_f, permutations=999, by = 'margin', strata = map_f[[strata_var]] )
    else
      permanova <- vegan::adonis( frm , data = map_f, permutations=999, by = 'margin')
    
    pos = which(prm_vars == var)
    
    pval[i] =  permanova$aov.tab$`Pr(>F)`[pos]
    R2[i] = permanova$aov.tab$R2[pos]
  }
  res = check_var_NA(vars_2_check_here, map)
  var_num = res[[1]]; NA_num = res[[2]]; samp_num = res[[3]]
  
  perm_df = data.frame(var = vars_2_check_here, R2 = R2, pval = pval, category_num = var_num, NA_num = NA_num, sample_number = samp_num)
  perm_df = perm_df[rev(order(perm_df$R2)),]
  
  perm_type = sprintf('ControlFor_%s',paste(extra_vars, collapse = '_'))
  if (!is.na(strata_var))
    perm_type = sprintf('%s_Strata_%s', perm_type, strata_var)
  if (length(extra_vars) == 0)
    perm_type = gsub('ControlFor_','',perm_type)
  perm_df$Permanova_extra_prms = perm_type
  
  return(perm_df)
}

## organize the data and than run permanova
permanova_organize_and_run = function( map, dist, vars_2_check, min_notNA_num = 3, extra_vars = c(), strata_var = NA )
{
  res = permanova_organize_data(map, dist,  vars_2_check, min_notNA_num)
  map = res[[1]]
  dist = res[[2]]
  vars_2_check_here = res[[3]]
  
  perm_df = permanova_run( map, dist, vars_2_check_here, extra_vars, strata_var  )
  return(perm_df)
}

## gets a list of file paths and merges them into one datafrmae
merge_df_files_to_one_df = function(fls, suffix = 'txt')
{
  fls = fls[grepl(suffix, fls)]
  na_str = c('no_data','_','NA','unknown', 'other','na','No_followup','ND','Unkown')
  map = read.table(fls[1], header = T, sep = '\t', na.strings = na_str)
  map$file = fls[1]
  if(length(fls) > 1)
  {
    for ( i in 2:length(fls))
    {
      temp = read.table(fls[i], header = T, sep = '\t', na.strings = na_str)
      temp$file = fls[i]
      map = rbind(map, temp)
    }
  }
  return(map)
}

## organize the permanova results data frame - 
## remove parameters with under filter_num samples or more than filter_num categories
## and set p value cutoffs and names
## also filtered only to paramters that were significant ( <0.1 ) in at least min_sig_to_plot settings
permanova_organize_res_df = function(df, vars_2_keep, filter_num = 10, min_sig_to_plot = 0, fdr_flag = F, filter_many_cats_flag = T)
{
  t = df
  t = t[t$var %in% vars_2_keep,]
  t$R2 = as.numeric(t$R2)
  
  ## filter mostly categoric (dates, etc.)
  t = t[t$sample_number>filter_num,]
  if (filter_many_cats_flag )
    t = t[ as.numeric(t$category_num)<filter_num | is.na(as.numeric(t$category_num)),]
  t = t[ t$category_num != t$sample_number, ]
  
  t$P.val = vector(mode = 'character',length = dim(t)[1])
  t$P.val[t$pval <= 0.1] = '<= 0.1'
  t$P.val[t$pval <= 0.05] = '<= 0.05'
  t$P.val[t$pval > 0.1] = '> 0.1'
  
  t$Type2 = sprintf('%s (n=%s)', t$Type, t$NA_num + t$sample_number)
  
  if ( min_sig_to_plot > 0 )
  {
   good_vars = c()
   for ( v in t$var )
   {
     pos = t$var == v
     if ( sum( t$P.val[pos] != '> 0.1') >= min_sig_to_plot )
     {
       good_vars = c(v, good_vars)
     }
   }
   good_vars = unique(good_vars)
   t = t[t$var %in% good_vars,]

   name = sprintf('%s_2sig', name)
  }
  return(t)
}

## organize the permanova results data frame - 
## remove parameters with under filter_num samples or more than filter_num categories
## and set p value cutoffs and names
## also filtered only to paramters that were significant ( <0.1 ) in at least min_sig_to_plot settings
permanova_organize_res_df_fdr = function(df, vars_2_keep, filter_num = 10, min_sig_to_plot = 0)
{
  t = df
  t = t[t$var %in% vars_2_keep,]
  t$R2 = as.numeric(t$R2)
  
  ## filter mostly categoric (dates, etc.)
  t = t[t$sample_number>filter_num,]
  t = t[ as.numeric(t$category_num)<filter_num | is.na(as.numeric(t$category_num)),]
  
  t$P.val = vector(mode = 'character',length = dim(t)[1])
  t$P.val[t$pval <= 0.1] = '<= 0.1'
  t$P.val[t$pval <= 0.05] = '<= 0.05'
  t$P.val[t$pval > 0.1] = '> 0.1'
  
  t$Type2 = sprintf('%s (n=%s)', t$Type, t$NA_num + t$sample_number)
  
  t = permanova_add_fdr(t)
  
  if ( min_sig_to_plot > 0 )
  {
    good_vars = c()
    for ( v in t$var )
    {
      pos = t$var == v
      if ( sum( t$Q.val[pos] != '> 0.25') >= min_sig_to_plot )
      {
        good_vars = c(v, good_vars)
      }
    }
    good_vars = unique(good_vars)
    t = t[t$var %in% good_vars,]
    
    name = sprintf('%s_2sig', name)
  }
  return(t)
}

permanova_add_fdr = function(t)
{
  t$qval = vector(mode = 'numeric',length = length(t$pval))
  for ( i in unique(t$Type2))
  {
    pos = t$Type2 == i
    t$qval[pos] = p.adjust(t$pval[pos], method = 'BH', n = sum(pos))
  }
  t$Q.val = vector(mode = 'character',length = dim(t)[1])
  t$Q.val[t$qval <= 0.25] = '<= 0.25'
  t$Q.val[t$qval <= 0.1] = '<= 0.1'
  t$Q.val[t$qval > 0.25] = '> 0.25'
  return(t)
}
