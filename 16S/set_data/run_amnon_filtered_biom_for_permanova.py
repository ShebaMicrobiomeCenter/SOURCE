#!/usr/bin/env python
import os
import sys
sys.path.insert(0, '/pita/users/tzipi/code/qiime/')
from qiime2_funcs import *

# source activate qiime2-2019.4

# setting veriables
verbose = True
trim_length = 150
threads = 20
decompres_flag = True
paired_flag = False

	
## merge files

# set variables
# map_file = 'data/SOURCE_16S_israel_china_map_v4.txt'
# out_path = 'res_israel_china_all_amnonFlt/'
# rar_num = 4000

# map_file = 'data/SOURCE_16S_israel_map_v4_Stool.txt'
# out_path = 'res_israel_stool_amnonFlt/'
# rar_num = 4000

# map_file = 'data/SOURCE_16S_israel_map_v4_noRS.txt'
# out_path = 'res_israel_noRS_amnonFlt/'
# rar_num = 4000

# map_file = 'data/SOURCE_16S_israel_map_v4.txt'
# out_path = 'res_israel_all_amnonFlt/'
# rar_num = 4000

# map_file = 'data/SOURCE_16S_china_map_v4.txt'
# out_path = 'res_china_all_amnonFlt//'
# rar_num = 33000

# map_file = 'data/SOURCE_16S_china_map_v4_Stool.txt'
# out_path = 'res_china_stool_amnonFlt//'
# rar_num = 33000

# map_file = 'data/SOURCE_16S_china_map_v4_Stool_rural.txt'
# out_path = 'res_china_stool_rural_amnonFlt//'
# rar_num = 33000

map_file = 'data/SOURCE_16S_china_map_v4_Stool_rural_ageMatched.txt'
out_path = 'res_china_stool_rural_ageMatched_amnonFlt//'
rar_num = 33000

# map_file = 'data/SOURCE_16S_china_map_v4_Stool_urban_and_CD.txt'
# out_path = 'res_china_stool_urban_and_CD_amnonFlt//'
# rar_num = 33000

# map_file = 'data/SOURCE_16S_china_map_v4_Stool_urban_and_CD_ageMatched.txt'
# out_path = 'res_china_stool_urban_and_CD_ageMatched_amnonFlt/'
# rar_num = 33000

os.system('mkdir ' + out_path)

# amnon_biom_file = 'data/SOURCE_16S_amnonFlt_israel_china_preprocessed.biom'
amnon_biom_file = 'data/SOURCE_israel_china_all_clean_amnon_mar31.biom'
os.system('qiime tools import \
	--input-path ' + amnon_biom_file + ' \
  	--type \'FeatureTable[Frequency]\' \
  	--input-format BIOMV210Format \
  	--output-path ' + out_path + '/table.qza')
# os.system('cp  ' + out_path + '/table.qza')
os.system('cp res_israel_china_all/rep-seqs.qza ' + out_path + '/rep-seqs.qza')

'''
# merge table
tables = ['/pita/users/tzipi/projects/16S/all_merges/DB1-19_merge/res/16S_DB1-19_merged/table.qza','/pita/users/tzipi/projects/16S/SOURCE_China/res/table.qza']
os.system('qiime feature-table merge \
	--i-tables ' + tables[0] + '\
	--i-tables ' + tables[1] + '\
	--o-merged-table ' + out_path + '/table.qza')

# merge rep-seqs
seqs = [ '/pita/users/tzipi/projects/16S/all_merges/DB1-19_merge/res/16S_DB1-19_merged/rep-seqs.qza','/pita/users/tzipi/projects/16S/SOURCE_China/res/rep-seqs.qza']
os.system('qiime feature-table merge-seqs \
	--i-data ' + seqs[0] + '\
	--i-data ' + seqs[1] + '\
	--o-merged-data ' + out_path + '/rep-seqs.qza')
'''

## filter samples by map, to remove samples not in the map (usually low barcode mistake samples here, because I used to the full barcodes maps)
os.system('cp ' + out_path + '/table.qza ' + out_path + '/table_all_samps.qza')
os.system('qiime feature-table filter-samples \
  --i-table ' + out_path + '/table_all_samps.qza \
  --m-metadata-file ' + map_file + ' \
  --o-filtered-table ' + out_path + '/table.qza')

## filtering general repseq to wanted samples
os.system('cp ' + out_path + '/rep-seqs.qza ' + out_path + '/rep-seqs_all_samps.qza')
os.system('qiime feature-table filter-seqs \
  --i-data ' + out_path + 'rep-seqs_all_samps.qza \
  --i-table ' + out_path + '/table.qza \
  --o-filtered-data ' + out_path + '/rep-seqs.qza')

## run downstream analysis
os.system('cp ' + map_file + ' ' + out_path + '/map.txt')

get_final_reads_count(out_path)
update_map_with_reads(out_path)
taxa_analysis(out_path, map_file, verbose, threads = threads)
# taxa_analysis_silva(out_path, map_file, verbose, threads = threads)
# phylogenetic_analysis(out_path, verbose, threads = threads)
phylogenetic_analysis_sepp(out_path, verbose, threads = threads)
diversity_analysis(out_path, map_file, rar_num, verbose, decompres_flag, threads = threads)

# clean_big_files(out_path)
open_files_for_R(out_path)

