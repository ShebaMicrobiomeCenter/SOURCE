library(ggplot2)

path = 'res/full_pairs/'

# cmp = 'Israel_metabolomicsMpink_Stool_vs_Israel_FFQ'
# xprm = 'Sucrose'
# # xprm = 'Histamine'
# # xprm = 'Lactose'
# # yprm = 'added_sugar_g_day'
# yprm = 'Fat_gr_day'
# # yprm = 'tryphtophan_mg_day'
# ylab = yprm
# # xprm = 'Hydroquinone'
# # yprm = 'Fish_and_poultry_servings_week'
# # ylab = 'Fish and poultry\nservings week'


# cmp = 'Israel_metabolomicsMblack_Stool_vs_Israel_metagenomics_species'
# xprm = 'N.Acetylasparagine'
# yprm = 'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Blautia|s__Blautia_faecis'
# ylab = 's__Blautia_faecis'
# xprm = 'N.Acetylasparagine'
# yprm = 'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Lacrimispora|s__Lacrimispora_amygdalina'
# ylab = 's__Lacrimispora_amygdalina'
# xprm = 'Melatonin'
# yprm = 'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Ruminococcaceae_unclassified|s__Ruminococcaceae_bacterium'
# ylab = 's__Ruminococcaceae_bacterium'



# xprm = 'Tryptophan'
# yprm = 'k__Bacteria|p__Firmicutes|c__Negativicutes|o__Veillonellales|f__Veillonellaceae|g__Veillonella|s__Veillonella_parvula'
# ylab = 's__Veillonella_parvula'
# xprm = 'Dopamine'
# yprm = 'k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Rikenellaceae|g__Alistipes|s__Alistipes_dispar'
# ylab = 's__Alistipes_dispar'
# xprm = 'Methionine'
# yprm = 'k__Bacteria|p__Firmicutes|c__Negativicutes|o__Veillonellales|f__Veillonellaceae|g__Veillonella|s__Veillonella_dispar'
# ylab = 's__Veillonella_dispar'

# cmp = 'Israel_metabolomicsMbrown_Stool_vs_Israel_metagenomics_species'
# xprm = 'Tyramine'
# yprm = 'k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Escherichia|s__Escherichia_coli'
# ylab = 's__Escherichia_coli'
# xprm = 'Malonate'
# yprm = 'k__Bacteria|p__Firmicutes|c__Negativicutes|o__Veillonellales|f__Veillonellaceae|g__Veillonella|s__Veillonella_parvula'
# ylab = 's__Veillonella_parvula'
# yprm = 'k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Pasteurellales|f__Pasteurellaceae|g__Haemophilus|s__Haemophilus_parainfluenzae'
# ylab = 's__Haemophilus_parainfluenzae'
# xprm = 'X3.Hydroxymethylglutarate'
# yprm = 'k__Bacteria|p__Firmicutes|c__Negativicutes|o__Veillonellales|f__Veillonellaceae|g__Veillonella|s__Veillonella_parvula'
# ylab = 's__Veillonella_parvula'
# xprm = 'Ureidopropionate'
# yprm = 'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Mediterraneibacter|s__Ruminococcus_gnavus'
# ylab = 's__Ruminococcus_gnavus'


# cmp = 'Israel_metabolomicsMbrown_Stool_vs_Israel_metagenomics_pathabundance'
# # xprm = 'Melatonin'
# # yprm = 'PWY-5130: 2-oxobutanoate degradation I'
# xprm = 'Ureidopropionate'
# yprm = 'BIOTIN-BIOSYNTHESIS-PWY: biotin biosynthesis I'

# cmp = 'Israel_metabolomicsMgreen_Stool_vs_Israel_metagenomics_pathabundance'
# xprm = 'Thiamine.Monophosphate'
# yprm = 'PWY-7953: UDP-N-acetylmuramoyl-pentapeptide biosynthesis III (meso-diaminopimelate containing)'
# ylab = 'PWY-7953'
# yprm = 'PWY-5686: UMP biosynthesis I'
# ylab = 'PWY-5686:\nUMP biosynthesis I'
# yprm = 'PWY-7790: UMP biosynthesis II'
# ylab = 'PWY-7790:\nUMP biosynthesis II'
# yprm = 'PWY-3841: folate transformations II (plants)'
# ylab = 'PWY-3841:\nfolate transformations II (plants)'

# cmp = 'Israel_metabolomicsMsalmon_Stool_vs_Israel_metagenomics_pathabundance'
# xprm = 'GABA'
# yprm = 'PWY-5104: L-isoleucine biosynthesis IV'
# ylab = 'PWY-5104: L-isoleucine\nbiosynthesis IV'

cmp = 'Israel_metabolomicsMtan_Stool_vs_Israel_metagenomics_ecsNamed'
xprm = 'cis.7.10.13.16.Docosatetraenoic.acid'
# xprm = 'X5.Z..8.Z..11.Z..Eicosatrienoic.Acid'
# yprm = 'NAD(+) diphosphatase (3.6.1.22)'
# ylab = 'NAD(+) diphosphatase (3.6.1.22)'
# yprm = 'Adenylate cyclase (4.6.1.1)'
# ylab = 'Adenylate cyclase (4.6.1.1)'
yprm = 'Enoyl-CoA hydratase (4.2.1.17)'
ylab = 'Enoyl-CoA hydratase (4.2.1.17)'

x_file = sprintf('%s/%s/x.txt', path, cmp)
y_file = sprintf('%s/%s/y.txt', path, cmp)

x=read.table(x_file, sep = '\t',comment.char = '', quote = '')
y=read.table(y_file, sep = '\t',comment.char = '', quote = '')

x=as.data.frame(t(x))
y=as.data.frame(t(y))

map = read.table(file = '../../metadata/alona_data/SOURCE_Israel_China_data_v11.txt', sep ='\t', header = T, row.names = 1)
map = map[row.names(x),]
map$Dx = gsub('healthy','Control',map$Dx)

p = ggplot(x, aes( x= x[[xprm]], y= y[[yprm]] )) +
  geom_point(size=3, aes(colour = map$Dx)) + theme_classic() +
  scale_colour_manual(name= 'Dx', values = c('#4D96FF','#FF6B6B'))+
  geom_smooth( method = 'lm', colour = 'gray', se=F) +
  scale_y_log10() +
  scale_x_log10() +
  xlab(xprm) + ylab(ylab)
p

temp = data.frame(x=x[[xprm]], y=y[[yprm]] )
lm_fit <- lm(y ~ x, data=temp)
lm_fit

reg<-lm(formula = y ~ x,
        data=temp)
coeff<-coefficients(reg)
intercept<-coeff[1]
slope<- coeff[2]

x_epsilon = min(temp$x[temp$x!=0])/10
temp$x[temp$x==0] = x_epsilon 
y_epsilon = min(temp$y[temp$y!=0])/10
temp$y[temp$y==0] = y_epsilon 

p = ggplot(temp, aes( x= x, y= y)) +
  geom_point(size=3, aes(colour = map$Dx)) + theme_classic() +
  scale_colour_manual(name= 'Dx', values = c('#4D96FF','#FF6B6B')) +
  geom_smooth(method = 'lm', colour = 'gray', se = F) + 
  # geom_abline(slope=lm_fit$coefficients[2],
  #             yintercept=lm_fit$coefficients[1],
  #             color="gray",
  #             linewidth=1) +
  # geom_line(data=data.frame(x=temp$x,y=intercept+slope*temp$x), color = 'gray',size=1) +
  scale_y_log10() +
  scale_x_log10() +
  xlab(xprm) + ylab(ylab)
p

out_file = sprintf('res/examples/%s_%s_%s.jpeg', cmp, xprm, yprm)
ggsave(filename = out_file, plot = p, device = 'jpeg',width = 4,height = 2.8)


names(y)[grepl('merd',names(y))]

