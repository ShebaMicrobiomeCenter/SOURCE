

PT=$'res_china_stool_rural_ageMatched_amnonFlt/'

# for beta diversity Patient_group2
qiime diversity beta-group-significance \
	--i-distance-matrix $PT/core-metrics-results/unweighted_unifrac_distance_matrix.qza \
	--m-metadata-file $PT/map.txt \
	--m-metadata-column Patient_group2 \
	--p-pairwise \
	--o-visualization $PT/core-metrics-results/unweighted-unifrac-Patient_group2-significance

qiime tools export \
	--input-path $PT/core-metrics-results/unweighted-unifrac-Patient_group2-significance.qzv \
	--output-path $PT/core-metrics-results/unweighted-unifrac-Patient_group2-significance'

