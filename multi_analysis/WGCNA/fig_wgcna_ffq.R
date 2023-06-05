## run the code to create the basic WGCNA FFQ heatmaps for israel and china
source("rnaSeq_main/WGCNA_using_other_modules_fdr.R")
out_path = 'rnaSeq_main/final_figs/'

col2 = colorRampPalette(rev( c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                               "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061") ))(200)
metadata_prms = c('Gender','Age_years','BMI',
  'Active_Smoker_manual','Previous_or_Current_Tobacco_Use_manual',
  'Disease','inflammed','CRP_numeric','Calprotectin_numeric',# 'CRP_log','Calprotectin_log',
  'health_index_same_sample','health_index_stool')

# set datafrmae for merged data
df_israel = res_israel[[2]]
df_china = res_china[[2]]

df_israel$Country = 'Israel'
df_china$Country = 'China'

df = rbind(df_israel, df_china)

df$var_type = ifelse(df$variable %in% metadata_prms, 'Metadata','Diet')

df$Country = factor(df$Country, levels = c('Israel','China'))
df$var_type = factor(df$var_type, levels = c('Metadata','Diet'))

df$id = sprintf('%s_%s',df$Country, df$variable)
dfm = as.data.frame(maditr::dcast(data = df[df$var_type =='Diet',],formula = id~module, 
                                  fun.aggregate = sum,value.var = "correlation") )
row.names(dfm) = dfm$id; dfm = dfm[,-1]
ht_israel_ffq = heatmap(as.matrix(dfm[grepl('Israel',row.names(dfm)),]), scale='none', col=col2, symm = F)
ht_China_ffq = heatmap(as.matrix(dfm[grepl('China',row.names(dfm)),]), scale='none', col=col2, symm = F)

df$id[df$var_type == 'Metadata'] = as.character(df$variable[df$var_type == 'Metadata'])
df$id[df$var_type == 'Metadata'] = as.character(df$variable[df$var_type == 'Metadata'])
df$id = factor(df$id, levels = c( rev(metadata_prms), 
                                  row.names(dfm[grepl('Israel',row.names(dfm)),])[ht_israel_ffq$rowInd],
                                  row.names(dfm[grepl('China',row.names(dfm)),])[ht_China_ffq$rowInd]
                                  ) )

## create the figure
library(ggplot2)
heatmap = ggplot(df, aes(x=module,y=id,fill = correlation)) +
  geom_tile() + xlab('Module') + 
  geom_tile(data = df[df$q_value <= min_q & !df$variable %in% metadata_prms,], 
            color = 'black', size = 1, aes(x=module,y=id,fill = correlation, alpha = q_value <= min_q)) +
  scale_fill_gradientn(colours = col2, limits = c(-1,1)) +
  scale_alpha_manual(values = c(0.2), name = sprintf('Q <= %s', min_q))+ 
  # geom_point(aes(shape = p_value<=0.05)) + scale_shape_manual(values = c(NA, 8)) +
  # geom_text(aes(label=sprintf('%.2f\n%.2e',correlation,p_value )), size=3) + 
  scale_colour_manual(values = c('black','white')) + 
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size= 10)) + 
  # theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size= 10)) + coord_flip() +
  scale_x_discrete(expand=c(0,0)) + 
  scale_y_discrete(expand=c(0,0)) + 
  scale_y_discrete(expand=c(0,0), 
                   labels = c('health_index_stool' = 'health_index_stool',
                              'health_index_same_sample' = 'health_index_same_sample',
                              'Calprotectin_numeric' = 'Calprotectin',
                              'CRP_numeric' = 'CRP',
                              'inflammed' = 'inflammation',
                              'Disease' = 'Disease',
                              'Previous_or_Current_Tobacco_Use_manual' = 'Previous_or_current_tobacco_use',
                              'Active_Smoker_manual' = 'Active_Smoking',
                              'BMI' = 'BMI',
                              'Age_years' = 'Age',
                              'Gender' = 'Gender',
                              'Israel_water_g_day' = 'Water',
                              'Israel_caffeine_mg_day' = 'Caffeine',
                              'Israel_X15i_ii._Additional_sugar_in._Coffee_teaspoons.lumps_cup__n' = 'Additional_sugar_in_coffee (qst)',
                              'Israel_Red_and_processed_meat_servings_week' = 'Red_and_processed_meat',
                              'Israel_Vitamin_B12_mcg_day' = 'Vitamin_B12',
                              'Israel_vitamin_k_mcg_day' = 'Vitamin_k',
                              'Israel_folate_mcg_day' = 'Folate',
                              'Israel_X15b._Vegetables_all_types__n' = 'Vegetables',
                              'Israel_vitamin_b6_mg_day' = 'Vitamin_b6',
                              'Israel_Protein_percentage_Kcal' = 'Protein_percentage_Kcal',
                              'Israel_Fish_and_poultry_servings_week' = 'Fish_and_poultry',
                              'Israel_Vitamin_D_mcg_day' = 'Vitamin_D',
                              'Israel_maganese_mg_day' = 'Maganese',
                              'Israel_omega_3_g_day' = 'Omega_3',
                              'Israel_Fat_percentage_Kcal' = 'Fat_percentage_Kcal',
                              'Israel_Olive_oil_servings_day' = 'Olive_oil',
                              'Israel_Fat_gr_day' = 'Fat',
                              'Israel_MUFA_g_day' = 'MUFA',
                              'Israel_Carbohydrates_percentage_Kcal'=  'Carbohydrates_percentage_Kcal',
                              'Israel_Carbohydrates_gr_day' = 'Carbohydrates',
                              'Israel_alpha_fat_acids_g_day' = 'Alpha_fat_acids',
                              'Israel_starch_g_day' = 'Starch',
                              'Israel_selenium_mcg_day' = 'Selenium',
                              'Israel_iodine_mcg_day' = 'Iodine',
                              'Israel_Olives.nuts.seed_servings_day' ='Olives_nuts_and_seeds',
                              'Israel_sugar_g_day' = 'Sugar',
                              'Israel_Dairy_servings_day' = 'Dairy',
                              'Israel_X15i_i._Additional_sugar_in.__choice.Breakfast_cereals_with_milk__n' = 'Additional_sugar_in_breakfast_cereals_with_milk',
                              'Israel_total_cis_g_day' = 'Total_cis',
                              'Israel_Saturated_fat_gr_day' = 'Saturated_fat',
                              'Israel_Saturated_fat_percentage_Kcal' = 'Saturated_fat_percentage_Kcal',
                              'Israel_Eggs_servings_week' = 'Eggs',
                              'Israel_X15d._Bread_slices_day__n' = 'Bread_slices (qst)',
                              'Israel_X15k._Drinks_._Juice__n' = 'Juice_drinking (qst)',
                              'Israel_X15l._Drinks_._Soft_drinks__n' ='Soft_drinks_drinking (qst)',
                              'Israel_added_sugar_g_day' = 'Added_sugar',
                              'Israel_Added_suger_percentage_Kcal' = 'Added_suger_percentage_Kcal',
                              'China_X15a._Fruit_all_type__n' = 'Fruit (qst)',
                              'China_Vitamin_D_mcg_day' ='Vitamin_D',
                              'China_X15b._Vegetables_all_types__n' = 'Vegetables (qst)',
                              'China_X15g._Breakfast_cereals_._Cornflakes.type__n' = 'Breakfast_cereals_cornflakes (qst)',
                              'China_Added_suger_percentage_Kcal' = 'Added_suger_percentage_Kcal',
                              'China_Fat_gr_day' = 'Fat',
                              'China_Saturated_fat_gr_day' = 'Saturated_fat',
                              'China_PUFA_g_day' = 'PUFA',
                              'China_X15j._Fast_food_e.g._food_from_a_hot.dog_stand_or_a_hamburger_restaurant__n'= 'Fast_food_e.g_hot_dog_or_hamburger (qst)')) + 
  ylab('') + guides(color = 'none', alpha) + 
  geom_text(aes(label=signif(p_value, 1), colour =p_value < 0.001 ), size=2.3) + 
  facet_grid(var_type + Country~., scales = 'free', space = 'free')

ggsave(sprintf('%s/module_cor_heatmap_israel_china_FFQ_p%s_q%s.pdf',out_path , min_p, min_q),
       plot = heatmap, device = 'pdf', width =9,height = 9)

write.table(x = df, 
            file = sprintf('%s/module_cor_heatmap_israel_china_FFQ_p%s_q%s.txt',out_path , min_p, min_q),
            quote = F, sep = '\t', row.names = F)
              

