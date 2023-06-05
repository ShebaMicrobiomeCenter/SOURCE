#!/usr/bin/env python

import os
import glob
import csv
import numpy
import pandas as pd
import numpy as np
import re

# merges fastq files (for big files or pulling of few samples)
# merges first list1 file with first list2 file under first id, etc.
def merge_fastqs(in_path, out_path, list1, list2, ids):
    # input check
    if len(list1) != len(list2) or len(list2) != len(ids):
        print('Input length does not fit')
        return
    for i in range(0,len(ids)):
        os.system('cat ' + in_path + list1[i] + ' ' + in_path + list2[i] + ' > ' + out_path + ids[i])       

# change name (and possibly path) of files
# save_flag = True to save original files, False to remove then
def change_names(in_names, out_names ,in_path = './', out_path = './', save_flag = True, add_txt = ''):
    if len(in_names) != len(out_names):
        print('Input length does not fit')
        return
    for i in range(0,len(in_names)):
        if (save_flag):
            os.system('cp ' + in_path + in_names[i] + ' ' + out_path + add_txt + out_names[i])
        else:
            os.system('mv ' + in_path + in_names[i] + ' ' + out_path + add_txt + out_names[i])

# run kallisto for a single fastq file by path
def kallisto_for_file_singleEnd(in_path, out_path):
    kallisto_index = '/pita/users/tzipi/DB/kallisto_index/hg_gencode_v24'
    #getting name
    np_file = in_path.split('/')[-1]
    file_name = np_file[0:-6]
    # kallisto
    os.system('kallisto quant -l 160 -s 20 -o ' +  out_path + file_name + ' -i ' + kallisto_index + ' --bias --single ' + in_path)
            
# have to make sure the files in in_path are compatiable = 
# all files have a match, and names fit untill ending _r1/_r2. 
# fits files by order, so if there is a problem the r1 r2 will not fit
def run_kallisto_paired_kal_42_5_gencode_v24(in_path, out_path, add_name = '', r1_suffix = '_r1.fastq', r2_suffix = '_r2.fastq', kallisto_index = '/pita/users/tzipi/DB/kallisto_index/hg_gencode_v24'):
    os.system('mkdir ' + out_path)

    # getting fastq files in in_folfer
    # files_r1 = glob.glob(in_path + '*_r1.fastq')
    files_r1 = glob.glob(in_path + '*' + r1_suffix)
    # files_r2 = glob.glob(in_path + '*_r2.fastq')
    files_r2 = glob.glob(in_path + '*' + r2_suffix)
    files_r1.sort()
    files_r2.sort()
    
    for i in range(0,len(files_r1)):
        #getting name
        np_file = files_r1[i].split('/')[-1]
        file_name1 = np_file[0:-len(r1_suffix)]
        
        #checking it is the same name as in r2
        np_file = files_r2[i].split('/')[-1]
        file_name2 = np_file[0:-len(r2_suffix)]
        if file_name1!=file_name2:
            print(file_name1 + ' != ' + file_name2)
        else:
            file_name1 = add_name + file_name1
            # kallisto
            os.system('kallisto quant -o ' +  out_path + file_name1 + ' -i ' + kallisto_index + ' --bias ' + files_r1[i] + ' ' + files_r2[i])
    
        print(file_name1)

def run_kallisto_paired_kal_42_5_gencode_v24_v2(in_path, out_path, add_name = '', r1_suffix = '_r1.fastq', r2_suffix = '_r2.fastq', kallisto_index = '/pita/users/tzipi/DB/kallisto_index/hg_gencode_v24', remove_name_start_flag = True):
    os.system('mkdir ' + out_path)

    # getting fastq files in in_folfer
    # files_r1 = glob.glob(in_path + '*_r1.fastq')
    files_r1 = glob.glob(in_path + '*' + r1_suffix)
    # files_r2 = glob.glob(in_path + '*_r2.fastq')
    files_r2 = glob.glob(in_path + '*' + r2_suffix)
    files_r1.sort()
    files_r2.sort()
    
    for i in range(0,len(files_r1)):

        #getting name
        np_file = files_r1[i].split('/')[-1]
        file_name1 = np_file[0:-len(r1_suffix)]
        
        #checking it is the same name as in r2
        np_file = files_r2[i].split('/')[-1]
        file_name2 = np_file[0:-len(r2_suffix)]
        if file_name1!=file_name2:
            print(file_name1 + ' != ' + file_name2)
        else:
            file_name1 = add_name + file_name1
            if remove_name_start_flag:
                file_name1 = re.sub(".*=", "", file_name1)
            # kallisto
            os.system('kallisto quant -o ' +  out_path + file_name1 + ' -i ' + kallisto_index + ' --bias ' + files_r1[i] + ' ' + files_r2[i])
    
        print(file_name1)


# running kallisto on single end fastq in folder. if there are files 
# to ignore they can be filtered by file name substrings in bad_files_txt
def run_kallisto_single_kal_42_5_gencode_v24(in_path, out_path, bad_files_txt = [], kallisto_index = '/pita/users/tzipi/DB/kallisto_index/hg_gencode_v24', s = 20):

    os.system('mkdir ' + out_path)

    # getting fastq files in in_folfer
    files = glob.glob(in_path + '*.fastq*')
    files.sort()
    
    # remove files containing bad text (usually paired end suffix)
    for bad in bad_files_txt:
        files = [k for k in files if bad not in k ]
    
    for file in files:
        #getting name
        np_file = file.split('/')[-1]
        file_name = np_file[0:-6]
        
        # kallisto
        os.system('kallisto quant -l 160 -s ' + str(s) + ' -o ' +  out_path + file_name + ' -i ' + kallisto_index + ' --bias --single ' + file)
        
        print(file_name)

# running kallisto on all files in the folder, single or paired end
# have to make sure the files in in_path are compatiable = 
# all files have a match, and names fit untill suffix. 
# fits files by order, so if there is a problem the r1 r2 will not fit
def run_kallisto_mixedEnd_kal_43_gencode_v24(in_path, out_path, add_name = '', r1_suffix = '_r1.fastq', r2_suffix = '_r2.fastq'):
    # runnign for paired end data
    run_kallisto_paired_kal_43_gencode_v24(in_path, out_path, add_name, r1_suffix, r2_suffix)
    # running for single end data, ignore paired end files by suffix
    suffixes = [r1_suffix, r2_suffix]
    run_kallisto_single_kal_43_gencode_v24(in_path, out_path, suffixes )
    
def sum_trans_2_gene_1file(in_file, out_file):
    
    id = [x[0] for x in csv.reader(open( in_file ,'r'),delimiter = '\t')]
    length = [x[1] for x in csv.reader(open( in_file ,'r'),delimiter = '\t')]
    eff_length = [x[2] for x in csv.reader(open( in_file ,'r'),delimiter = '\t')]
    est_counts = [x[3] for x in csv.reader(open( in_file ,'r'),delimiter = '\t')]
    tpm = [x[4] for x in csv.reader(open( in_file ,'r'),delimiter = '\t')]

    del id[0]
    del length[0]
    del eff_length[0]
    del est_counts[0]
    del tpm[0]

    gene_name = [x.split('|')[5] for x in id ]
    gene_type = [x.split('|')[7] for x in id ]
    
    # removing anything after an underscore ('_') from gene names
    for i in range(0,len(gene_name)):
        if '_' in gene_name[i]:
            gene_name[i] = gene_name[i].split('_')[0]
            
    u_gene  = list(set(gene_name))
    u_gene.sort()
    
    u_length = [None]*len(u_gene)
    u_eff_length = [None]*len(u_gene)
    u_est_counts = [None]*len(u_gene)
    tpm_sum = [None]*len(u_gene)
    u_type = [None]*len(u_gene)
    transcripts_count = [None]*len(u_gene)
    for i in range(0,len(u_gene)):
        u_length[i] = 0
        u_eff_length[i] = 0
        u_est_counts[i] = 0
        tpm_sum[i] = 0
        transcripts_count[i] = 0
        for j in range(0,len(gene_name)):
            if gene_name[j] == u_gene[i]:
                u_length[i] = max(u_length[i], float(length[j]) )
                u_eff_length[i] = max(u_eff_length[i], float(eff_length[j]) )
                u_est_counts[i] = u_est_counts[i] + float(est_counts[j])
                tpm_sum[i] = tpm_sum[i] + float(tpm[j])
                u_type[i] = gene_type[j]
                transcripts_count[i] = transcripts_count[i] + 1 

    of = open(out_file,'w')
    of.write('target_id\tlength\teff_length\test_counts\ttpm\ttype\ttrans_count\n')
    
    for i in range(0,len(u_gene)):
        of.write( u_gene[i] + '\t' + str(u_length[i]) + '\t' + str(u_eff_length[i]) +  '\t' + str(u_est_counts[i])  + '\t' + str(tpm_sum[i]) + '\t' + u_type[i] + '\t' + str(transcripts_count[i]) + '\n')

    of.close()
    

def sum_trans_2_gene_1file_v2(in_file, out_file):
    df = pd.read_csv(in_file, sep = '\t')

    gene_name = [x.split('|')[5] for x in df['target_id'] ]
    gene_type = [x.split('|')[7] for x in df['target_id'] ]
    
    # removing anything after an underscore ('_') from gene names
    for i in range(0,len(gene_name)):
        if '_' in gene_name[i]:
            gene_name[i] = gene_name[i].split('_')[0]
            
    u_gene  = list(set(gene_name))
    u_gene.sort()
    
    u_length = [None]*len(u_gene)
    u_eff_length = [None]*len(u_gene)
    u_est_counts = [None]*len(u_gene)
    tpm_sum = [None]*len(u_gene)
    u_type = [None]*len(u_gene)
    transcripts_count = [None]*len(u_gene)
    for i in range(0,len(u_gene)):
        u_length[i] = 0
        u_eff_length[i] = 0
        u_est_counts[i] = 0
        tpm_sum[i] = 0
        transcripts_count[i] = 0
        pos = np.where([u_gene[i] == s for s in gene_name])[0]
        for j in range(0,len(pos)):
            u_length[i] = max(u_length[i], float( df['length'][ pos[j] ]) )
            u_eff_length[i] = max(u_eff_length[i], float( df['eff_length'][ pos[j] ]) )
            u_est_counts[i] = u_est_counts[i] + float( df['est_counts'][ pos[j] ])
            tpm_sum[i] = tpm_sum[i] + float( df['tpm'][ pos[j] ])
            u_type[i] = gene_type[ pos[j] ]
            transcripts_count[i] = transcripts_count[i] + 1 

    of = open(out_file,'w')
    of.write('target_id\tlength\teff_length\test_counts\ttpm\ttype\ttrans_count\n')
    for i in range(0,len(u_gene)):
        of.write( u_gene[i] + '\t' + str(u_length[i]) + '\t' + str(u_eff_length[i]) +  '\t' + str(u_est_counts[i])  + '\t' + str(tpm_sum[i]) + '\t' + u_type[i] + '\t' + str(transcripts_count[i]) + '\n')

    of.close()

def trans2gene(in_path, out_path, rev_falg = False):
    os.system('mkdir ' + out_path)

    paths = glob.glob(in_path + '*/')
    paths.sort()
    # reverse order of files to run. good for long runs in parallal from 2 directions
    if rev_falg:
        paths.reverse()
        
    for path in paths:
        new_id = path.split('/')[-2]
        print(new_id)
        
        in_file = in_path + new_id + '/abundance.tsv'
        out_file = out_path + new_id + '_abundance.tsv'
        sum_trans_2_gene_1file(in_file, out_file)
        print(out_file)

# used to make files match naming sence that works with my other code. 
# highly recommended to backup data before use
def rename_fastq_read_mark_to_suffix(in_path, read_mark = 'read1_',rep_mark = '', old_suffix = '[.]fastq', new_suffix = '_r1.fastq'):
    files = glob.glob(in_path + '/*' + read_mark + '*')
    old_files = files
    import re
    files = [re.sub(read_mark,rep_mark, x) for x in files]
    files = [re.sub(old_suffix,new_suffix, x) for x in files]
    for i in range(len(files)):
        os.system('mv ' + old_files[i] + ' ' + files[i])


