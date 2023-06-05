DATA_PATH=$'data/'
OUT_PATH=$'qiime2_res/'

BIOM_NAME=$'SOURCE2022_Stool_Data_MZM13_N2_filtered_v2'
DIR_NAME=$'SOURCE2022_Stool_Data_MZM13_N2_filtered_v2'
MAP_FILE=$DATA_PATH'/SOURCE2022_Stool_metadata.txt'

BIOM_FILE=$DATA_PATH/$BIOM_NAME'.biom'
BIOM_TEXT_FILE=$DATA_PATH/$BIOM_NAME'.tsv'


mkdir $OUT_PATH
OUT_PATH=$OUT_PATH$DIR_NAME'/'
mkdir $OUT_PATH

cp $MAP_FILE $OUT_PATH'/map.tsv'

# biom convert -i $BIOM_FILE -o $BIOM_TEXT_FILE --to-tsv
# mkdir $DATA_PATH'/rep-seq/'
# cut -f1 $BIOM_TEXT_FILE | awk -F"\t" 'NR > 2 {print ">"$1"\n"$1}' > $DATA_PATH'/rep-seq/dna-sequences.fasta'

biom convert -i $BIOM_TEXT_FILE -o $BIOM_FILE --to-json

# importing rge biom amnon send into qiime2 -
# both biom and rep-seq (made the fansta manually
qiime tools import \
  --input-path $BIOM_FILE \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV100Format \
  --output-path $OUT_PATH'/'$BIOM_NAME'_biom_full.qza'
 

# qiime tools import \
# 	--type FeatureData[Sequence] \
# 	--input-format DNASequencesDirectoryFormat \
# 	--input-path $DATA_PATH/rep-seq/ \
# 	--output-path $DATA_PATH/rep-seq.qza
	
## filtering general biom to wanted samples
qiime feature-table filter-samples \
  --i-table $OUT_PATH'/'$BIOM_NAME'_biom_full.qza' \
  --m-metadata-file $MAP_FILE \
  --o-filtered-table $OUT_PATH'/'$BIOM_NAME'_biom.qza' 

# cp $DATA_PATH/rep-seq.qza $OUT_PATH/rep-seq.qza
		
# qiime fragment-insertion sepp \
# 	--i-representative-sequences $DATA_PATH/rep-seq.qza \
# 	--o-tree $OUT_PATH/rooted-tree.qza  \
# 	--o-placements $OUT_PATH/insertion-placements.qza \
# 	--p-threads 10

# qiime diversity core-metrics \
# 	--i-table $OUT_PATH/table.qza \
# 	--p-sampling-depth 9999 \
# 	--m-metadata-file $DATA_PATH/$NAME'_sample.txt' \
# 	--output-dir $OUT_PATH/core-metrics-results \
#	--p-n-jobs 20

mkdir $OUT_PATH/core-metrics-results

DISTACE_TYPE=$'canberra'
# DISTACE_TYPE=$'braycurtis'

qiime diversity beta \
	--i-table $OUT_PATH'/'$BIOM_NAME'_biom.qza' \
	--p-metric $DISTACE_TYPE \
	--p-n-jobs 20 \
	--o-distance-matrix $OUT_PATH'/core-metrics-results/'$BIOM_NAME'_'$DISTACE_TYPE'_distance_matrix'

qiime diversity pcoa \
	--i-distance-matrix $OUT_PATH'/core-metrics-results/'$BIOM_NAME'_'$DISTACE_TYPE'_distance_matrix.qza' \
	--o-pcoa $OUT_PATH'/core-metrics-results/'$BIOM_NAME'_'$DISTACE_TYPE'_pcoa_results'

qiime emperor plot \
	--i-pcoa $OUT_PATH'/core-metrics-results/'$BIOM_NAME'_'$DISTACE_TYPE'_pcoa_results.qza' \
	--m-metadata-file $MAP_FILE \
	--p-ignore-missing-samples \
	--o-visualization $OUT_PATH'/core-metrics-results/'$BIOM_NAME'_'$DISTACE_TYPE'_emperor'

qiime tools export \
 	--input-path $OUT_PATH'/core-metrics-results/'$BIOM_NAME'_'$DISTACE_TYPE'_emperor.qzv' \
 	--output-path $OUT_PATH'/core-metrics-results/'$BIOM_NAME'_'$DISTACE_TYPE'_emperor'
 
## make pcoa-biplot
qiime feature-table relative-frequency \
	--i-table $OUT_PATH'/'$BIOM_NAME'_biom.qza' \
	--o-relative-frequency-table $OUT_PATH'/'$BIOM_NAME'_table_RA'

qiime diversity pcoa-biplot \
	--i-pcoa $OUT_PATH'/core-metrics-results/'$BIOM_NAME'_'$DISTACE_TYPE'_pcoa_results.qza' \
	--i-features $OUT_PATH'/'$BIOM_NAME'_table_RA.qza' \
	--o-biplot $OUT_PATH'/core-metrics-results/'$BIOM_NAME'_'$DISTACE_TYPE'_pcoa_biplot_results' 

qiime emperor biplot \
	--i-biplot $OUT_PATH'/core-metrics-results/'$BIOM_NAME'_'$DISTACE_TYPE'_pcoa_biplot_results.qza' \
	--m-sample-metadata-file $MAP_FILE \
	--o-visualization $OUT_PATH'/core-metrics-results/'$BIOM_NAME'_'$DISTACE_TYPE'_pcoa_biplot_emperor' 

qiime tools export \
 	--input-path $OUT_PATH'/core-metrics-results/'$BIOM_NAME'_'$DISTACE_TYPE'_pcoa_biplot_emperor.qzv'  \
 	--output-path $OUT_PATH'/core-metrics-results/'$BIOM_NAME'_'$DISTACE_TYPE'_pcoa_biplot_emperor' 


qiime tools export \
 	--input-path $OUT_PATH'/core-metrics-results/'$BIOM_NAME'_'$DISTACE_TYPE'_distance_matrix.qza'  \
 	--output-path $OUT_PATH'/core-metrics-results/'$BIOM_NAME'_'$DISTACE_TYPE'_distance_matrix' 


