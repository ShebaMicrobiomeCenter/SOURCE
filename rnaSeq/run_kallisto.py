#!/usr/bin/env python
import os
import pandas as pd 
import sys
sys.path.insert(0, '/pita/users/tzipi/code/rnaSeq')
from run_kallisto_funcs import *

kallisto_index = '/pita/users/tzipi/bin/kallisto_linux-v0.42.5/index/hg_gencode_v24'

# setting variables (paths)
in_path = '/mnt/storage/pub/tzipi/SOURCE_China/organized/Transcriptome/include/'

out_path = 'kallisto/'
os.system('mkdir ' + out_path)

# run_kallisto_paired_kal_42_5_gencode_v24(in_path, out_path, r1_suffix = '_R1.fastq.gz', r2_suffix = '_R2.fastq.gz', kallisto_index = kallisto_index)

run_kallisto_paired_kal_42_5_gencode_v24(in_path, out_path, r1_suffix = '_1.fq.gz', r2_suffix = '_2.fq.gz', kallisto_index = kallisto_index)
'''
# merging transcript level to gene level
gene_out_path = out_path + 'kallisto_gene/'
os.system('mkdir ' + gene_out_path)
trans2gene(kal_out_path, gene_out_path)
'''

