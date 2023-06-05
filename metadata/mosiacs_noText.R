library(ggplot2)
library(ggmosaic)
library(patchwork)

ratio = 4
cols = c('#E57C23','#1A4D2E') 

## metadata
metadata_map_file = '/pita/users/tzipi/projects/multiomics/SOURCE/metadata/alona_data/SOURCE_Israel_China_data_v11.txt'
map = read.table(metadata_map_file, header = T, sep = '\t')

map$Patient_group3 = map$Patient_group2
map$Patient_group3[map$Patient_group3 == 'Chinese_crohns'] = 'China CD (n=40)'
map$Patient_group3[map$Patient_group3 == 'Urban_health'] = 'China Urban (n=121)'
map$Patient_group3[map$Patient_group3 == 'Rural_health_>50%_in_city'] = 'China Rural-Urban (n=74)'
map$Patient_group3[map$Patient_group3 == 'Rural_health_<50%_in_city'] = 'China Rural (n=88)'
map$Patient_group3[map$Patient_group3 == 'Israeli_Crohns'] = 'Israel CD (n=25)'
map$Patient_group3[map$Patient_group3 == 'Israeli_healthy'] = 'Israel Control (n=32)'

map$Patient_group4 = map$Patient_group3
map$Patient_group4[map$Patient_group3 %in% 
                     c('China Rural-Urban (n=74)',
                       'China Rural (n=88)')] = 'China Rural (n=162)'


prm = 'X17d._Flush_toilet_'
prm_name = 'Flush Toilet'
prm_lvls = c('Yes','No')

groups = c('China CD (n=40)','China Urban (n=121)')

map[[prm]] = factor(map[[prm]], levels =prm_lvls)

map_f = map[map$Patient_group3 %in% groups,]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(X17d._Flush_toilet_ , Patient_group3), 
                  fill = X17d._Flush_toilet_), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) + 
  scale_fill_manual(values = cols, name = '')

groups = c('China Rural-Urban (n=74)','China Rural (n=88)')
map_f = map[map$Patient_group3 %in% groups,]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g2 = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(X17d._Flush_toilet_, Patient_group3), 
                  fill = X17d._Flush_toilet_), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  scale_fill_manual(values = cols, name = '')

groups = c('Israel CD (n=25)','Israel Control (n=32)')
map_f = map[map$Patient_group3 %in% groups &!is.na(map$X17d._Flush_toilet_),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g3 = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(X17d._Flush_toilet_, Patient_group3), 
                  fill = X17d._Flush_toilet_), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  scale_fill_manual(values = cols, name = '')

groups = c('China Rural (n=162)','China Urban (n=121)')
map_f = map[map$Patient_group4 %in% groups &!is.na(map$X17d._Flush_toilet_),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g4 = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(X17d._Flush_toilet_, Patient_group4), 
                  fill = X17d._Flush_toilet_), colour = 'black')+
  theme_mosaic() + 
  xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  scale_fill_manual(values = cols, name = '')

m_flush = (g2 + theme(legend.position = 'none')) + 
  (g4 + theme(legend.position = 'none') + ylab('')) + 
  (g + theme(legend.position = 'none') + ylab(''))  +
  (g3 + ylab('')) + plot_layout(nrow = 1)

##### soft drink
prm_name = 'Soft drinks'
prm = 'Soft_drinks'
prm_lvls = c('Weekly or more','Less than weekly')
map$Soft_drinks = ifelse(map$X15l._Drinks_._Soft_drinks_ == 'Less_frequently','Less than weekly','Weekly or more')
# map$Soft_drinks = factor(map$Soft_drinks, levels = c('Weekly or more','Less than weekly'))
groups = c('China CD (n=40)','China Urban (n=121)')

map[[prm]] = factor(map[[prm]], levels =prm_lvls)

map_f = map[map$Patient_group3 %in% groups,]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(Soft_drinks , Patient_group3), fill = Soft_drinks), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) + 
  scale_fill_manual(values = cols, name = '')

groups = c('China Rural-Urban (n=74)','China Rural (n=88)')
map_f = map[map$Patient_group3 %in% groups,]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g2 = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(Soft_drinks, Patient_group3), 
                  fill = Soft_drinks), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  scale_fill_manual(values = cols, name = '')

groups = c('Israel CD (n=25)','Israel Control (n=32)')
map_f = map[map$Patient_group3 %in% groups &!is.na(map$Soft_drinks),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g3 = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(Soft_drinks, Patient_group3), 
                  fill = Soft_drinks), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  scale_fill_manual(values = cols, name = '')

groups = c('China Rural (n=162)','China Urban (n=121)')
map_f = map[map$Patient_group4 %in% groups &!is.na(map$Soft_drinks),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g4 = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(Soft_drinks, Patient_group4), 
                  fill = Soft_drinks), colour = 'black')+
  theme_mosaic() + 
  xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  scale_fill_manual(values = cols, name = '')

m_soft = (g2 + theme(legend.position = 'none')) + 
  (g4 + theme(legend.position = 'none') + ylab('')) + 
  (g + theme(legend.position = 'none') + ylab(''))  +
  (g3 + ylab('')) + plot_layout(nrow = 1)



##### Tap water

prm_name = 'Tap water'
prm = 'Tap_water'
prm_lvls = c('Yes','No')
map$Tap_water = map$X17a._In_house_water_tap_
groups = c('China CD (n=40)','China Urban (n=121)')

map[[prm]] = factor(map[[prm]], levels =prm_lvls)

map_f = map[map$Patient_group3 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(Tap_water , Patient_group3), 
                  fill = Tap_water), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) + 
  scale_fill_manual(values = cols, name = '')

groups = c('China Rural-Urban (n=74)','China Rural (n=88)')
map_f = map[map$Patient_group3 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g2 = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(Tap_water, Patient_group3), 
                  fill = Tap_water), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  scale_fill_manual(values = cols, name = '')

groups = c('Israel CD (n=25)','Israel Control (n=32)')
map_f = map[map$Patient_group3 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g3 = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(Tap_water, Patient_group3), 
                  fill = Tap_water), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  scale_fill_manual(values = cols, name = '')

groups = c('China Rural (n=162)','China Urban (n=121)')
map_f = map[map$Patient_group4 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g4 = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(Tap_water, Patient_group4), 
                  fill = Tap_water), colour = 'black')+
  theme_mosaic() + 
  xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  scale_fill_manual(values = cols, name = '')

m_tap = (g2 + theme(legend.position = 'none')) + 
  (g4 + theme(legend.position = 'none') + ylab('')) + 
  (g + theme(legend.position = 'none') + ylab(''))  +
  (g3 + ylab('')) + plot_layout(nrow = 1)


##### farm animals

prm_name = 'Farm animals'
prm = 'Farm_animals'
prm_lvls = c('No','Yes')
map$Farm_animals = map$X10g._Farm_animals
groups = c('China CD (n=40)','China Urban (n=121)')

map[[prm]] = factor(map[[prm]], levels =prm_lvls)

map_f = map[map$Patient_group3 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(Farm_animals , Patient_group3), 
                  fill = Farm_animals), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) + 
  scale_fill_manual(values = cols, name = '')

groups = c('China Rural-Urban (n=74)','China Rural (n=88)')
map_f = map[map$Patient_group3 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g2 = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(Farm_animals, Patient_group3), 
                  fill = Farm_animals), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  scale_fill_manual(values = cols, name = '')

groups = c('Israel CD (n=25)','Israel Control (n=32)')
map_f = map[map$Patient_group3 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g3 = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(Farm_animals, Patient_group3), 
                  fill = Farm_animals), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  scale_fill_manual(values = cols, name = '')

groups = c('China Rural (n=162)','China Urban (n=121)')
map_f = map[map$Patient_group4 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g4 = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(Farm_animals, Patient_group4), 
                  fill = Farm_animals), colour = 'black')+
  theme_mosaic() + 
  xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  scale_fill_manual(values = cols, name = '')

m_farm = (g2 + theme(legend.position = 'none')) + 
  (g4 + theme(legend.position = 'none') + ylab('')) + 
  (g + theme(legend.position = 'none') + ylab(''))  +
  (g3 + ylab('')) + plot_layout(nrow = 1)


##### coffee

prm_name = 'Coffee'
prm = 'Coffee'
prm_lvls = c('>1 cup/day','None')
map$Coffee = map$X15m._Drinks_._Coffee_
map$Coffee = ifelse(map$Coffee == '0','None','>1 cup/day')
groups = c('China CD (n=40)','China Urban (n=121)')

map[[prm]] = factor(map[[prm]], levels =prm_lvls)

map_f = map[map$Patient_group3 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(Coffee , Patient_group3), 
                  fill = Coffee), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) + 
  scale_fill_manual(values = cols, name = '')

groups = c('China Rural-Urban (n=74)','China Rural (n=88)')
map_f = map[map$Patient_group3 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g2 = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(Coffee, Patient_group3), 
                  fill = Coffee), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  scale_fill_manual(values = cols, name = '')

groups = c('Israel CD (n=25)','Israel Control (n=32)')
map_f = map[map$Patient_group3 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g3 = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(Coffee, Patient_group3), 
                  fill = Coffee), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  scale_fill_manual(values = cols, name = '')

groups = c('China Rural (n=162)','China Urban (n=121)')
map_f = map[map$Patient_group4 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g4 = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(Coffee, Patient_group4), 
                  fill = Coffee), colour = 'black')+
  theme_mosaic() + 
  xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  scale_fill_manual(values = cols, name = '')

m_coffee = (g2 + theme(legend.position = 'none')) + 
  (g4 + theme(legend.position = 'none') + ylab('')) + 
  (g + theme(legend.position = 'none') + ylab(''))  +
  (g3 + ylab('')) + plot_layout(nrow = 1)

##### rural 5

prm_name = 'Rural before age 5'
prm = 'rural5'
prm_lvls = c('No','Yes')
map$rural5 = map$X16a._Infant_0.5years_
map$rural5 = ifelse(map$rural5 %in% c('Rural','Countryside'),'Yes','No')
groups = c('China CD (n=40)','China Urban (n=121)')

map[[prm]] = factor(map[[prm]], levels =prm_lvls)

map_f = map[map$Patient_group3 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(rural5 , Patient_group3), 
                  fill = rural5), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) + 
  scale_fill_manual(values = cols, name = '')

groups = c('China Rural-Urban (n=74)','China Rural (n=88)')
map_f = map[map$Patient_group3 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g2 = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(rural5, Patient_group3), 
                  fill = rural5), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  scale_fill_manual(values = cols, name = '')

groups = c('Israel CD (n=25)','Israel Control (n=32)')
map_f = map[map$Patient_group3 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g3 = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(rural5, Patient_group3), 
                  fill = rural5), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  scale_fill_manual(values = cols, name = '')

groups = c('China Rural (n=162)','China Urban (n=121)')
map_f = map[map$Patient_group4 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g4 = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(rural5, Patient_group4), 
                  fill = rural5), colour = 'black')+
  theme_mosaic() + 
  xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  scale_fill_manual(values = cols, name = '')

m_rural5 = (g2 + theme(legend.position = 'none')) + 
  (g4 + theme(legend.position = 'none') + ylab('')) + 
  (g + theme(legend.position = 'none') + ylab(''))  +
  (g3 + ylab('')) + plot_layout(nrow = 1)


##### smoking f

prm_name = 'Current smoking (females)'
prm = 'smoke_f'
prm_lvls = c('Yes','No')
map$smoke_f = map$X12.Active_Smoker_manual
map_old = map
map = map[map$Gender=='female',]

groups = c('China CD (n=40)','China Urban (n=121)')

map[[prm]] = factor(map[[prm]], levels =prm_lvls)

map_f = map[map$Patient_group3 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(smoke_f , Patient_group3), 
                  fill = smoke_f), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) + 
  scale_fill_manual(values = cols, name = '')

groups = c('China Rural-Urban (n=74)','China Rural (n=88)')
map_f = map[map$Patient_group3 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g2 = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(smoke_f, Patient_group3), 
                  fill = smoke_f), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  scale_fill_manual(values = cols, name = '')

groups = c('Israel CD (n=25)','Israel Control (n=32)')
map_f = map[map$Patient_group3 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g3 = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(smoke_f, Patient_group3), 
                  fill = smoke_f), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  scale_fill_manual(values = cols, name = '')

groups = c('China Rural (n=162)','China Urban (n=121)')
map_f = map[map$Patient_group4 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g4 = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(smoke_f, Patient_group4), 
                  fill = smoke_f), colour = 'black')+
  theme_mosaic() + 
  xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  scale_fill_manual(values = cols, name = '')

m_smoke_f = (g2 + theme(legend.position = 'none')) + 
  (g4 + theme(legend.position = 'none') + ylab('')) + 
  (g + theme(legend.position = 'none') + ylab(''))  +
  (g3 + ylab('')) + plot_layout(nrow = 1)
map = map_old

##### smoking m

prm_name = 'Current smoking (males)'
prm = 'smoke_m'
prm_lvls = c('Yes','No')
map$smoke_m = map$X12.Active_Smoker_manual
map_old = map
map = map[map$Gender=='male',]

groups = c('China CD (n=40)','China Urban (n=121)')

map[[prm]] = factor(map[[prm]], levels =prm_lvls)

map_f = map[map$Patient_group3 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(smoke_m , Patient_group3), 
                  fill = smoke_m), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) + 
  scale_fill_manual(values = cols, name = '')

groups = c('China Rural-Urban (n=74)','China Rural (n=88)')
map_f = map[map$Patient_group3 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g2 = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(smoke_m, Patient_group3), 
                  fill = smoke_m), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  scale_fill_manual(values = cols, name = '')

groups = c('Israel CD (n=25)','Israel Control (n=32)')
map_f = map[map$Patient_group3 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g3 = ggplot(data = map_f) +
  geom_mosaic(aes(x = product(smoke_m, Patient_group3), 
                  fill = smoke_m), colour = 'black')+
  theme_mosaic() + xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  scale_fill_manual(values = cols, name = '')

groups = c('China Rural (n=162)','China Urban (n=121)')
map_f = map[map$Patient_group4 %in% groups & !is.na(map[[prm]]),]
prm_vals = factor(map_f[[prm]], levels =prm_lvls)
g4 = ggplot(data = map_f) + 
  geom_mosaic(aes(x = product(smoke_m, Patient_group4), 
                  fill = smoke_m), colour = 'black')+
  theme_mosaic() + 
  xlab('') + ylab(prm_name)  +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=15),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio=ratio) +
  # axis.text.x = element_text(angle =45, hjust = 0), aspect.ratio=ratio) +
  # scale_x_productlist(position = "top") +
  scale_fill_manual(values = cols, name = '')

m_smoke_m = (g2 + theme(legend.position = 'none')) + 
  (g4 + theme(legend.position = 'none') + ylab('')) + 
  (g + theme(legend.position = 'none') + ylab(''))  +
  (g3 + ylab('')) + plot_layout(nrow = 1)
map = map_old


library(ggpubr)

g = ggarrange(m_flush, m_tap, m_rural5, m_farm, m_soft, m_coffee, m_smoke_m, m_smoke_f,
              ncol = 2, nrow = 4, align = 'hv')

ggsave( 'mosiacs_noText.pdf',plot = g, device = 'pdf', width =12,height = 12)




