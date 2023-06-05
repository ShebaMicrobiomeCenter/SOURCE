#!/usr/bin/env python
import os
import sys
sys.path.insert(0, '/pita/users/tzipi/code/qiime/')
from qiime2_funcs import *

# source activate qiime2-2019.1

# setting veriables
verbose = True
trim_length = 150
threads = 10
rar_num = 33000
decompres_flag = True
paired_flag = False


## setting parameteres
# input_path = '/mnt/storage/pub/tzipi/SOURCE_China/organized/16S/fastq_split/'


# map_file = '/mnt/storage/pub/tzipi/SOURCE_China/organized/16S/metadata/SOURCE_china_16S_premap.txt'
# manifest_file = '/mnt/storage/pub/tzipi/SOURCE_China/organized/16S/metadata/SOURCE_china_16S_SE_manifest.txt'

input_path = '/mnt/storage/pub/tzipi/SOURCE_China/organized/16S/fastq_split_v4/'
out_path = 'res_stool/'

map_file = 'data/SOURCE_china_16S_stool_premap.txt'
manifest_file = 'data/SOURCE_china_16S_SE_stool_v4_manifest.txt'

os.system('mkdir ' + out_path)
os.system('cp ' + map_file + ' ' + out_path + '/map.txt')

## running qiime basic analysis
# import_data_to_qiime2( input_path, out_path, verbose, paired_flag = paired_flag )
# demultiplex( out_path, map_file, verbose, paired_flag = paired_flag )

## manifest import (no need to demultiplex here)
'''
manifest_file = 'data/manifest.txt'
os.system('qiime tools import \
              --type \'SampleData[SequencesWithQuality]\' \
              --input-path ' + manifest_file + ' \
              --output-path ' + out_path + '/demux.qza \
              --input-format PairedEndFastqManifestPhred33')
'''

import_data_with_manifest(input_path, out_path, verbose, paired_flag, manifest_file)

# DADA2( out_path, 0, trim_length, verbose, threads = threads, paired_flag = paired_flag )
deblur( out_path, map_file, trim_length, verbose, threads = threads )	
# phylogenetic_analysis(out_path, verbose, threads = threads)
phylogenetic_analysis_sepp(out_path, verbose, threads = threads)
diversity_analysis(out_path, map_file, rar_num, verbose, decompres_flag, threads = threads)
taxa_analysis(out_path, map_file, verbose, threads = threads)
get_final_reads_count(out_path)
# clean_big_files(out_path)
update_map_with_reads(out_path)
