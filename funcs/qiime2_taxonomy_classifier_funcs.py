# creates and Greenegens classifier specificaly fitted to our data - 
# the 515F/806R primmer set and 175bp long reads.
# the Greenegene version is 8_13 and identity is 99
def GG_classifier_train(classifier):
	# load dowloaded GG_13_8 data
	os.system('qiime tools import \
  --type \'FeatureData[Sequence]\' \
  --input-path 99_otus.fasta \
  --output-path gg_otus.qza')

	os.system('qiime tools import \
  --type \'FeatureData[Taxonomy]\' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path 99_otu_taxonomy.txt \
  --output-path gg_ref-taxonomy.qza')

	# fitting the classifier to our specific data - 
	# Cutting to read length and fitting primers. 
	# The primers used in the qiime2 moving pictures tutorial are 515F/806R, the same as we used here.
	os.system('qiime feature-classifier extract-reads \
  --i-sequences gg_otus.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-trunc-len 150 \
  --o-reads gg_ref-seqs.qza')
	
	# training the classifier
	os.system('qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads gg_ref-seqs.qza \
  --i-reference-taxonomy gg_ref-taxonomy.qza \
  --o-classifier ' + classifier)
	
	
	
# creates and silva classifier specificaly fitted to our data - 
# the 515F/806R primmer set and 150bp long reads.
# using silva release 138 qza file from the qiime2 resource page: 
# https://docs.qiime2.org/2020.11/data-resources/
def GG_classifier_train_silva(classifier):
	# fitting the classifier to our specific data - 
	# Cutting to read length and fitting primers. 
	# The primers used in the qiime2 moving pictures tutorial are 515F/806R, the same as we used here.
	os.system('qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-trunc-len 150 \
  --o-reads silva-138-99-seqs_150bp.qza')
	
	os.system('qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs-515-806.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-trunc-len 150 \
  --o-reads silva-138-99-seqs-515-806_150bp.qza')
	
	# training the classifier
	os.system('qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-99-seqs_150bp.qza \
  --i-reference-taxonomy silva-138-99-tax.qza \
  --o-classifier ' + classifier)
	
