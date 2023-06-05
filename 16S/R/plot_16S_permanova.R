source('/pita/users/tzipi/code/R_figs/permanova_funcs.R')

# set parameters
min_sig_to_plot = 2
filter_num = 10
sig2_flag = T

path = 'res/16S_un_all'
out_path = 'res_final/'

## read data
fls = list.files(path = path)
fls = sprintf('%s/%s',path, fls)
df = merge_df_files_to_one_df(fls)

# set data
df$Type = gsub('.*/','',df$file)
df$Type = gsub('[.]txt','',df$Type)

bad_vars = c('Chronic_comorbidities','Concomitant_medications',
             'X19a._Death_of_family_member._','Calprotectin_ug_g',
             'Version_of_ICF_signed','Sibling_2','kcal_per_kweight')
vars_2_keep =  sort(unique(df$var))
vars_2_keep =vars_2_keep[! vars_2_keep %in% bad_vars]

# write results to table for supp.
## need to add smoking
df_text = permanova_organize_res_df(df, vars_2_keep, filter_num = 6, min_sig_to_plot = 0)
write.table(x = df_text, file = sprintf('%s/16S_un_all_permanova_table.txt', out_path), quote = F, sep = '\t',row.names = F)

t = permanova_organize_res_df(df, vars_2_keep, filter_num = 6, min_sig_to_plot = min_sig_to_plot)
t$Type2

ffq_file = '../../metadata/FFQ/FFQ_v10.txt'
ffq = read.table(ffq_file, header = T, sep = '\t', row.names = 1)
wanted_vars = c('Age_years','age','Gender','BMI','Weight_kilograms',
                'Height_cm','pn_ID','Smoking',
                'Regular_physical_activities_walking_jogging_cycling_swimming_over_30_minutes_or_similar_activities_at_present',
                'Regular_physical_activities_walking_jogging_cycling_swimming_over_30_minutes_or_similar_activities_before_diagnosis_of_IBD',
                'Dx','Patient_group','Patient_group2','Disease_Status',
                'inflammation_status',
                'Calprotectin_numeric','Calprotectin_ug_g',
                'CRP_numeric','CRP_Absolute','CRP_Category','CRP_mg_L',
                'FCP_Category','FCP_Absolute','CRP_Category',
                'Faith_pd','Dysbiosis_Index','Health_Index',
                'health_index_same_sample','health_index_stool')
# t = t[t$var %in% c(wanted_vars, names(ffq)), ]
# wanted_vars = c('X4._Do_you_have_or_have_had_long_lasting_repetitive_problems_with_your_stomach.',
#                 'X9c._Have_you_received_any_antibiotics_in_the_last_3_months.',
#                 'X1._Do_you_have_Siblings','X6i._Asthma_')
# t = t[t$var %in% c(wanted_vars), ]
t$var_old = t$var

t$var[t$var == 'X10e._Aquarium_fishes_'] = 'Aquarium fish'
t$var[t$var == 'X8a._Measles'] = 'Past measles infection'
t$var[t$var == 'X15k._Drinks_._Juice__n'] = 'Juice drinking'
t$var[t$var == 'X15l._Drinks_._Soft_drinks__n'] = 'Soft drinks drinking'
t$var[t$var == 'X6i._Asthma_'] = 'Asthma'
t$var[t$var == 'Age_years'] = 'Age'
t$var[t$var == 'X4._Do_you_have_or_have_had_long_lasting_repetitive_problems_with_your_stomach.'] = 'Long lasting_stomach problems'
t$var[t$var == 'X9c._Have_you_received_any_antibiotics_in_the_last_3_months.'] = 'Antibiotics last 3 months'
t$var[t$var == 'X10g._Farm_animals'] = 'Farm animals'
t$var[t$var == 'X1._Do_you_have_Siblings'] = 'Having siblings'

t$var = gsub('_g_day','',t$var)
t$var = gsub('_mg_day','',t$var)
t$var = gsub('_gr_day','',t$var)
t$var = gsub('_servings_day','',t$var)
t$var = gsub('_servings_week','',t$var)
t$var = gsub('Olives.nuts.seed','Olives nuts and seeds',t$var)
t$var = gsub('_',' ',t$var)

if (min_sig_to_plot == 2 & path == 'res/16S_un_all')
{
  prm_ord = c('Added suger percentage Kcal','sugar',
              'Carbohydrates percentage Kcal','Fruits',
              'Fibers',
              'Dairy','Calcium','Sodium',
              'Legumes',
              'Olives nuts and seeds','Fat',
              'Saturated fat','PUFA','MUFA',
              'Ferrus','Protein','Protein percentage Kcal',
              'Red and processed meat',
              'Energy Kcal day',# 'kcal_per_kweight',
              'Juice drinking','Soft drinks drinking',
              'Age','BMI',
              'Asthma','Past measles infection',
              'Long lasting stomach problems','Antibiotics last 3 months',
              'Aquarium fish','Farm animals','Having siblings')
  temp  = unique(t$var)
  t$var = factor(t$var, levels = rev( c(prm_ord, temp[!temp %in% prm_ord]) ) )
}
 
t$var_type = 'Questionnaire'
t$var_type[t$var_old %in% names(ffq) & (!grepl('X15',t$var_old) )] = 'FFQ'
t$Country = ifelse(grepl('China',t$Type),'China','Israel')

t$group = 'Control'
t$group[ grepl('CD',t$Type) ] = 'CD'
t$group[ grepl('Urban',t$Type) ] = 'Urban'
t$group[ grepl('Rural>50%',t$Type) ] = 'Rural-Urban'
t$group[ grepl('Rural<50%',t$Type) ] = 'Rural'
t$Type3 = sprintf('%s %s %s', t$Country, t$group, gsub('.* ','',t$Type2))
t$Type3 = factor(t$Type3, levels =c('China CD (n=40)','China Urban (n=121)',
                                    'China Rural-Urban (n=74)','China Rural (n=88)',
                                    'Israel CD (n=20)','Israel Control (n=28)'))
# plot
library(viridis)

g = ggplot(t, aes(y=var, x=Type3, fill = R2*100)) + 
  geom_tile(colour = 'black') + 
  # geom_point(shape = 8, aes(colour = P.val))+
  geom_point(colour = 'white' , aes(shape = P.val, alpha = P.val)) +
  scale_shape_manual(values=c(8, 4, 46), breaks = c('<= 0.05','<= 0.1','> 0.1')) +
  scale_alpha_manual(values = c(1,1,0), breaks = c('<= 0.05','<= 0.1','> 0.1')) +#  guides(alpha=FALSE) + 
  scale_colour_manual(values = c('black',NA), na.translate=FALSE) + ## NA TO GRAY
  theme_classic() + ylab('') + xlab('') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  # scale_fill_viridis() + 
  # scale_fill_viridis(option = 'D',
  scale_fill_gradient(low = 'lightblue',high='blue4') +
  # scale_fill_gradient(low = 'lightblue',high='blue4',
  #   rescaler = function(x, to = c(0, 1), from = NULL) 
  # { ifelse(x<6,  scales::rescale(x,to = to, from = c(min(x, na.rm = TRUE), 10)), 1) } ) +
  scale_x_discrete(expand=c(0,0)) +
  theme(legend.key = element_rect(fill = "blue4")) +
  labs(fill = 'Variance %\nexplained', colour = 'P.val <= 0.1')+  
  scale_y_discrete(expand=c(0,0)) +
  facet_grid(var_type~Country, scales = 'free',space = 'free') + 
  theme(panel.background = element_rect(fill = 'gray', colour = 'gray'),
        panel.spacing = unit(0,'lines'))

h = 2+length(unique(t$var))/3
w = 3+length(unique(t$Type))
# # # ggsave(sprintf('%s/permanova_fig_minSig%s.pdf', path, min_sig_to_plot),plot = g, device = 'pdf', width = w,height = h, limitsize = FALSE)
# # ggsave(sprintf('%s/permanova_fig_minSig%s.pdf', path, min_sig_to_plot),plot = g, device = 'pdf', width =5,height = 6.5, limitsize = FALSE)
# # ggsave(sprintf('%s/permanova_fig_minSig%s_env.pdf', path, min_sig_to_plot),plot = g, device = 'pdf', width =8 ,height = 9, limitsize = FALSE)
# ggsave(sprintf('%s/permanova_fig_minSig%s_env.pdf', out_path, min_sig_to_plot),plot = g, device = 'pdf', width =w ,height = h, limitsize = FALSE)
ggsave(sprintf('%s/permanova_fig_minSig%s_env.pdf', out_path, min_sig_to_plot),
       plot = g, device = 'pdf', width =5 ,height =5.5, limitsize = FALSE)

g2 = ggplot(t, aes(y=var, x=Type, fill = R2*100)) +
  geom_tile(colour = 'black') +
  # geom_point(shape = 8, aes(colour = P.val))+
  # geom_point(colour = 'white' , aes(shape = P.val, alpha = P.val)) +
  # scale_shape_manual(values=c(8, 4, 46), breaks = c('<= 0.05','<= 0.1','> 0.1')) +
  # scale_alpha_manual(values = c(1,1,0), breaks = c('<= 0.05','<= 0.1','> 0.1')) +#  guides(alpha=FALSE) +
  # scale_colour_manual(values = c('black',NA), na.translate=FALSE) + ## NA TO GRAY
  theme_classic() + ylab('') + xlab('') +
  # geom_text(aes(label=sprintf('%.2f%%\n%.3f',R2*100,pval ), color = R2*100>5), size=2)+
  geom_text(aes(label=sprintf('%.2f%%\n%.3f',R2*100,pval ), color = pval<0.05), size=2)+
  scale_color_manual(values = c('black','white')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # scale_fill_viridis() +
  # scale_fill_viridis(option = 'D',
  scale_fill_gradient(low = 'lightblue',high='blue4') +
  # scale_fill_gradient(low = 'lightblue',high='blue4',
  #   rescaler = function(x, to = c(0, 1), from = NULL)
  # { ifelse(x<6,  scales::rescale(x,to = to, from = c(min(x, na.rm = TRUE), 10)), 1) } ) +
  scale_x_discrete(expand=c(0,0)) +
  theme(legend.key = element_rect(fill = "blue4")) +
  # labs(fill = 'Variance %\nexplained', colour = 'P.val <= 0.1')+
  labs(fill = 'Variance %\nexplained')+
  scale_y_discrete(expand=c(0,0)) +
  theme(panel.background = element_rect(fill = 'gray', colour = 'gray'))

g2

# ggsave(sprintf('%s/permanova_fig_nums_minSig%s.pdf', path, min_sig_to_plot),plot = g2, device = 'pdf', width = w,height = h, limitsize = FALSE)

# 
