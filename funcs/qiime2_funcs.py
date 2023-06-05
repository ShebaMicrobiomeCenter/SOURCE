#!/usr/bin/env python
import os
import glob
import pandas as pd

# source activate qiime2-2019.4

## built for qiime2-2019.4. may noy work for other verisions

## updated to 2021.4

# source activate deblurenv
# deblur workflow --seqs-fp seqs.fna --output-dir deblur -w -t 150 -a 4 &
# source deactivate

# open the qzv files of unweighted_unifrac_emperor.qzv and faith-pd-group-significance.qzv
def abDiv_decompress(path):
	os.system('qiime tools export \
	--input-path ' + path + 'unweighted_unifrac_emperor.qzv \
	--output-path ' + path + 'unweighted_unifrac_emperor')
	os.system('qiime tools export \
	--input-path ' + path + 'faith-pd-group-significance.qzv \
	--output-path ' + path + 'faith-pd-group-significance')

# opens all qzv (vizualisation) files in the directory
def decompress_qzv(path):
	files = glob.glob(path + '/*.qzv')
	for f in files :
		name = f.split('/')[-1].split('.')[0]
		os.system('qiime tools export \
		--input-path ' + f + ' \
		--output-path ' + path + '/' + name)

def decompress_by_path(file_path):
	name = file_path.split('.')[0]
	os.system('qiime tools export \
	--input-path ' + file_path + '\
	--output-path ' + name )


# removing the big files after a qiime2 run to avoid cluttering the system
def clean_big_files(path):
	os.system('rm ' + path + '/emp-single-end-sequences.qza')
	os.system('rm ' + path + '/demux.qza')
	os.system('rm ' + path + '/demux-filtered.qza')
	
# import the fastq files into a qza file qiime2 can work with
# the forward.fastq.gz, reverse.fastq.gz and barcodes.fastq.gz files need to be in the 
# input_path folder (qiime2 needs the files to have those specific names)
def import_data_to_qiime2(input_path, out_path, verbose = False, paired_flag = False):
	if (paired_flag):	
		os.system('qiime tools import \
  		--type EMPPairedEndSequences \
  		--input-path ' + input_path + ' \
  		--output-path ' + out_path + '/emp-paired-end-sequences.qza')
	else:
		os.system('qiime tools import \
  		--type EMPSingleEndSequences \
  		--input-path ' + input_path + ' \
  		--output-path ' + out_path + '/emp-single-end-sequences.qza')
	if verbose:
		print('import complete')
	return

# demultiplex and pair if paired-end reads
def demultiplex(out_path, map_file, verbose = False, paired_flag = False ):
	# demultiplex
	if (paired_flag):
		os.system('qiime demux emp-paired \
		--i-seqs ' + out_path + '/emp-paired-end-sequences.qza \
		--m-barcodes-file ' + map_file + ' \
		--p-rev-comp-mapping-barcodes \
		--m-barcodes-column BarcodeSequence \
		--o-per-sample-sequences ' + out_path + 'demux.qza \
		--o-error-correction-details ' + out_path + 'demux-details.qza')
	else:
		os.system('qiime demux emp-single \
		--i-seqs ' + out_path + '/emp-single-end-sequences.qza \
		--m-barcodes-file ' + map_file + ' \
		--p-rev-comp-mapping-barcodes \
		--m-barcodes-column BarcodeSequence \
		--o-per-sample-sequences ' + out_path + 'demux.qza \
		--o-error-correction-details ' + out_path + 'demux-details.qza')
	# creates a .qzv visualization file of the results.
	# the visualizations can be viewed using the qiime2 
	# "qiime tools view" command
	os.system('qiime demux summarize \
 		--i-data ' + out_path + '/demux.qza \
 		--o-visualization ' + out_path + '/demux.qzv')	  
	# os.system('qiime tools view ' + out_path + 'demux.qzv')
	if verbose:
		print('demultiplex complete')
	return

# creating a manifest file for all the fastq.gz files in a given directory. 
def create_manifest_singleEnd(input_path, manifest_file):
	pwd = os.getcwd()
	full_path = pwd + '/' + input_path
	files = glob.glob(full_path + '/*.fastq.gz')
	ids = [fl.split('/')[-1].split('.')[0] for fl in files]
	temp={'sample-id': ids, 'absolute-filepath': files, 'direction': 'forward'}
	df = pd.DataFrame(data = temp)
	df.to_csv(manifest_file, index=False)

# if fastq files are already demultiplexed (usually downloaded data)
# doing the import using manifest. creates a manifest file (currently works for single end only)
# and uses it to create a demux.qza file. 
def import_data_with_manifest(input_path, out_path, verbose, paired_flag = False, manifest_file = ''):
	if (paired_flag):
		if (manifest_file == ''): 
			print('Does not support automatic paired end manifest file import yet.')
			return
		os.system('qiime tools import \
		--type \'SampleData[SequencesWithQuality]\' \
		--input-path ' + manifest_file + ' \
		--output-path ' + out_path + '/demux.qza \
		--input-format PairedEndFastqManifestPhred33V2')
	else:
		if (manifest_file == ''): 
			manifest_file = out_path + '/manifest.txt' 
			create_manifest_singleEnd(input_path, manifest_file)
		os.system('qiime tools import \
            --type \'SampleData[SequencesWithQuality]\' \
            --input-path ' + manifest_file + ' \
            --output-path ' + out_path + '/demux.qza \
            --input-format SingleEndFastqManifestPhred33V2')
	# creates a .qzv visualization file of the results.
	# the visualizations can be viewed using the qiime2 
	# "qiime tools view" command
	os.system('qiime demux summarize \
    --i-data ' + out_path + '/demux.qza \
    --o-visualization ' + out_path + '/demux.qzv')	 
	if verbose:
		print('import manifest demutliplexed data complete')


# cleaning data before using deblur
def deblur_clean(out_path):
	# initial quality filtering	
	'''
	os.system('qiime quality-filter q-score \
	--i-demux ' + out_path + '/demux.qza \
	--o-filtered-sequences ' + out_path + '/demux-filtered.qza \
	--o-filter-stats ' + out_path + '/demux-filter-stats.qza')

	# more stringent, like old qiime1 cleaning (this removes ~30% of reads)
	os.system('qiime quality-filter q-score \
	--i-demux ' + out_path + '/demux.qza \
	--o-filtered-sequences ' + out_path + '/demux-filtered.qza \
	--o-filter-stats ' + out_path + '/demux-filter-stats.qza \
	--p-min-quality 19 \
	--p-quality-window 1 \
	--p-min-length-fraction 0.75')
	'''
	# more stringent, less than old qiime1 cleaning 
	os.system('qiime quality-filter q-score \
	--i-demux ' + out_path + '/demux.qza \
	--o-filtered-sequences ' + out_path + '/demux-filtered.qza \
	--o-filter-stats ' + out_path + '/demux-filter-stats.qza \
	--p-min-quality 19 \
	--p-quality-window 3 \
	--p-min-length-fraction 0.10') # no need to remove here, deblur will remove it if shorter than the wanted trim_length.

# run deblur - do intiall cleaning, deblur and generate summary statistics
def deblur(out_path, map_file, trim_length = 150, verbose = False, threads = 10 ):
	deblur_clean(out_path)
	# run deblur
	os.system('qiime deblur denoise-16S \
	--i-demultiplexed-seqs ' + out_path + '/demux-filtered.qza \
	--p-trim-length ' + str(trim_length) + ' \
	--o-representative-sequences ' + out_path + '/rep-seqs.qza \
	--o-table ' + out_path + 'table.qza \
	--p-sample-stats \
	--o-stats ' + out_path + 'deblur-stats.qza \
	--p-no-hashed-feature-ids \
	--p-jobs-to-start ' + str(threads) )
	# generate summary statistics
	os.system('qiime metadata tabulate \
	--m-input-file ' + out_path + '/demux-filter-stats.qza \
	--o-visualization ' + out_path + '/demux-filter-stats.qzv')
	os.system('qiime deblur visualize-stats \
	--i-deblur-stats ' + out_path + '/deblur-stats.qza \
	--o-visualization ' + out_path + '/deblur-stats.qzv')
	if verbose:
		print('deblur complete')
	return

# runs the demultiplexed data through DADA2.
# This performs quality control, trunculates reads and removes bad
# quality reads and chimeras. Creates a feature table with clean SVs
# and a list of representive sequences (for all SVs that passed quality control)
def DADA2(out_path, trim_left = 0, trunc_len = 150, verbose = False, paired_flag = False, threads = 10):
	# Runs DADA2 - quality control + creates SVs table and list
	if (paired_flag):	
		os.system('qiime dada2 denoise-paired \
		--i-demultiplexed-seqs ' + out_path  + '/demux.qza \
		--p-trim-left-f ' + str(trim_left) + ' \
		--p-trunc-len-f ' + str(trunc_len) + ' \
		--p-trim-left-r ' + str(trim_left) + ' \
		--p-trunc-len-r ' + str(trunc_len) + ' \
		--o-representative-sequences ' + out_path + '/rep-seqs.qza \
		--o-table ' + out_path + '/table.qza \
		--o-denoising-stats ' + out_path + '/stats-dada2.qza \
		--p-no-hashed-feature-ids \
		--p-n-threads ' + str(threads) )
	else:
		os.system('qiime dada2 denoise-single \
		--i-demultiplexed-seqs ' + out_path  + '/demux.qza \
		--p-trim-left ' + str(trim_left) + ' \
		--p-trunc-len ' + str(trunc_len) + ' \
		--o-representative-sequences ' + out_path + '/rep-seqs.qza \
		--o-table ' + out_path + '/table.qza \
		--o-denoising-stats ' + out_path + '/stats-dada2.qza \
		--p-no-hashed-feature-ids \
		--p-n-threads ' + str(threads) )
	# creates a .qzv visualization files of the results.
	os.system('qiime metadata tabulate \
	--m-input-file ' + out_path + '/stats-dada2.qza \
	--o-visualization ' + out_path + '/stats-dada2.qzv')
	if verbose:
		print('DADA2 complete')
	return


# performes a phylogenetic analysis of the SVs. The resulting tree will be used 
# to calculate core-metrics of various alpha and beta diversity measures.
# each sample is rarefied to rar_num reads at this stage to avoid sample size bias
def phylogenetic_analysis(out_path, verbose = False, threads  = 10):
	# phylogenetic diversity analysis
	os.system('qiime phylogeny align-to-tree-mafft-fasttree \
	--i-sequences ' + out_path + '/rep-seqs.qza \
	--o-alignment ' + out_path + '/aligned-rep-seqs.qza \
	--o-masked-alignment ' + out_path + '/masked-aligned-rep-seqs.qza \
	--o-tree ' + out_path + '/unrooted-tree.qza \
	--o-rooted-tree ' + out_path + '/rooted-tree.qza \
	--p-n-threads ' + str(threads) )
	if verbose:
		print('phylogenetic analysis complete')
	return


def phylogenetic_analysis_sepp(out_path, verbose = False, threads = 10):
	# phylogenetic diversity analysis
	# run sepp on the available sequences
	## get database as reccomended bu qiime2 guide: (dowload to /pita/users/tzipi/DB/greengenes
	# wget \
	# -O "sepp-refs-gg-13-8.qza" \
	#   "https://data.qiime2.org/2019.10/common/sepp-refs-gg-13-8.qza"
	os.system('qiime fragment-insertion sepp \
		--i-representative-sequences ' + out_path + '/rep-seqs.qza \
		--i-reference-database /pita/users/tzipi/DB/greengenes/sepp-refs-gg-13-8.qza \
		--o-tree ' + out_path + '/rooted-tree.qza  \
		--o-placements ' + out_path + '/insertion-placements.qza \
		--p-threads ' + str(threads) )
	# removing sequences that were not inserted into the tree by sepp (weird sequnces not close to tree at all)
	os.system('qiime fragment-insertion filter-features \
		--i-table ' + out_path + '/table.qza \
		--i-tree ' + out_path + '/rooted-tree.qza \
		--o-filtered-table ' + out_path + '/table-filtered.qza \
		--o-removed-table ' + out_path + '/table-removed.qza \
		--verbose')
	# changing names so that the main talbe will be the filtered table
	os.system( 'mv ' + out_path + '/table.qza ' + out_path + '/table-unfiltered.qza' )
	os.system( 'mv ' + out_path + '/table-filtered.qza ' + out_path + '/table.qza' )
	if verbose:
		print('sepp complete')
	return

# performes alpha and beta diversity analysis
# each sample is rarefied to rar_num reads at this stage to avoid sample size bias
def diversity_analysis(out_path, map_file, rar_num = 2000, verbose = False, decompres_flag = True, threads = 10):
	# phylogenetic diversity analysis
	# if core-metrics-results folder already exists, move to old to avoid no overide problem
	# if os.path.isdir(out_path + 'core-metrics-results'):
	# 	os.system('mv ' + out_path + '/core-metrics-results ' + out_path + '/core-metrics-results_old')
	# perform basic diversity analysis
	os.system('qiime diversity core-metrics-phylogenetic \
	--i-phylogeny ' + out_path + '/rooted-tree.qza \
	--i-table ' + out_path + '/table.qza \
	--p-sampling-depth ' + str(rar_num) + ' \
	--m-metadata-file ' + map_file + ' \
	--output-dir ' + out_path + '/core-metrics-results \
	--p-n-jobs-or-threads ' + str(threads) )
	# check alpha diversity connection to variables
	os.system('qiime diversity alpha-group-significance \
	  --i-alpha-diversity ' + out_path + '/core-metrics-results/faith_pd_vector.qza \
	  --m-metadata-file ' + map_file + ' \
	  --o-visualization ' + out_path + '/core-metrics-results/faith-pd-group-significance.qzv')
	os.system('qiime diversity alpha-correlation \
	  --i-alpha-diversity ' + out_path + '/core-metrics-results/faith_pd_vector.qza \
	  --m-metadata-file ' + map_file + ' \
	  --o-visualization ' + out_path + '/core-metrics-results/faith-pd-correlation.qzv')
	if (decompres_flag):
		decompress_qzv(out_path + '/core-metrics-results/')
	if verbose:
		print('diversity analysis complete')
	return

# classifies SVs to taxonomy using the pre-trained classifier. makes taxonomy tables at all levels 
# and taxonomic interactive barplots (taxa-bar-plots.qzv)	
def taxa_analysis(out_path, map_file, verbose = False, classifier = '/pita/users/tzipi/DB/greengenes/gg_13_8_otus/13_8_i99_150bp_classifier/gg_99_classifier_150bp.qza', decompres_flag = True, threads = 10):
	# classifies feathures (SVs) to taxa using the pre-trained classifier
	os.system('qiime feature-classifier classify-sklearn \
	--i-classifier ' + classifier + ' \
	--i-reads ' + out_path + '/rep-seqs.qza \
	--o-classification ' + out_path + '/taxonomy.qza \
	--p-n-jobs ' + str(threads) )
	# creating a visualization of the resulting table (.qzv file)
	os.system('qiime metadata tabulate \
	--m-input-file ' + out_path + '/taxonomy.qza \
	--o-visualization ' + out_path + '/taxonomy.qzv')
	# making an interactive taxonomic barplot
	os.system('qiime taxa barplot \
	--i-table ' + out_path + '/table.qza \
	--i-taxonomy ' + out_path + '/taxonomy.qza \
	--m-metadata-file ' + map_file + ' \
	--o-visualization ' + out_path + '/taxa-bar-plots.qzv')
	if decompres_flag:
		decompress_qzv(out_path)
	if verbose:
		print('taxa analysis complete')
	return

# exports the biom table and creates a summary text file of the number
# of reads per sample. those reads should have passed all filteres.
def get_final_reads_count(out_path, decompres_flag = True):
	os.system('qiime tools export \
	--input-path ' + out_path + '/table.qza \
	--output-path ' + out_path + '/biom')
	os.system('biom summarize-table \
	-i ' + out_path + '/biom/feature-table.biom \
	-o ' + out_path + '/biom/biom_sum.txt')
	os.system('qiime feature-table summarize \
	--i-table ' + out_path + '/table.qza \
	--m-sample-metadata-file ' + out_path + '/map.txt \
	 --o-visualization ' + out_path + '/biom/table_sum.qzv')
	if decompres_flag:
		decompress_qzv(out_path + '/biom/')


# Does the upstream analysis from fastq files to biom
def qiime2_fastq_to_table(input_path, out_path, map_file, trim_length = 150, verbose = True, threads = 10 ):
	import_data_to_qiime2( input_path, out_path, verbose )
	demultiplex( out_path, map_file, verbose )
	# DADA2( out_path, 0, trim_length, verbose, threads = threads )
	deblur( out_path, map_file, trim_length, verbose, threads = threads )	

# Does the downstream analysis from table to alpha beta diversity results.
def qiime2_table_to_diversity_taxa(out_path, map_file, rar_num, verbose = True, decompres_flag =True, threads = 10):
	# phylogenetic_analysis(out_path, verbose, threads = threads)
	phylogenetic_analysis_sepp(out_path, verbose, threads = threads)
	diversity_analysis(out_path, map_file, rar_num, verbose, decompres_flag, threads = threads)
	taxa_analysis(out_path, map_file, verbose, threads = threads)

# add number of final reads per sample to map
def update_map_with_reads(out_path):
	import csv
	reads_file = out_path + '/biom/table_sum/sample-frequency-detail.csv'
	reads_df = pd.read_csv(reads_file, sep=',', names = ['SampleID','Reads_count'] )
	map_file = out_path + '/map.txt'
	map = pd.read_csv(map_file, sep='\t', header=0)
	map_reads = pd.merge(reads_df, map, how = 'right', on=['SampleID','SampleID'])
	with open(out_path + '/biom/map_reads.txt','w') as write_tsv:
	    write_tsv.write(map_reads.to_csv(sep='\t', na_rep='NA', index=False))

# decompressing the files needed for automatic R analysis.
def open_files_for_R(out_path):
	file_path = out_path + '/core-metrics-results/faith_pd_vector.qza'
	decompress_by_path(file_path)
	file_path = out_path + '/core-metrics-results/unweighted_unifrac_distance_matrix.qza'
	decompress_by_path(file_path)
	file_path = out_path + '/core-metrics-results/unweighted_unifrac_pcoa_results.qza'
	decompress_by_path(file_path)
