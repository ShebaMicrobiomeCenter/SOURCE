library(tximport)

# library('DESeq2')
source('/pita/users/tzipi/code/rnaSeq/DESeq2_DE_funcs.r')

two_groups = c('Crohns','healthy')
prm = 'Dx'

# location_prm = 'location_inflamation'
# # location_groups  = c('TI1','TI2','R3','R4')
# location_groups  = c('TI2','R4')
# # location_group = 'TI2'

location_prm = 'Location'
location_groups  = c('TI','R')

one_sample_per_patient_filter_flag = T

metadata_file = 'data/rnaSeq_SOURCE_china_map_v3.txt'

df = read.table(metadata_file, header = T, sep = '\t')

name = 'SOURCE_China'
input_dir = 'kallisto/'

# make a transcript to gene invertion table based on a random kallisto-output file
example_file = sprintf('%s/B202TI2/abundance.tsv', input_dir)
exp = read.table(example_file,sep="\t", header=T)
trans_names = exp$target_id

gene_name = as.character(trans_names)
ens_name = as.character(trans_names)
for ( i in 1:length(gene_name) )  {   gene_name[i] = strsplit(gene_name[i], '|', fixed = TRUE)[[1]][6] }
for ( i in 1:length(ens_name) )  {   ens_name[i] = strsplit(ens_name[i], '|', fixed = TRUE)[[1]][2] }
GENEID = sprintf('%s_%s',gene_name, ens_name)
# gene_name = gsub(pattern = '_.*',replacement = '',x = gene_name)
tx2gene = data.frame(TXNAME = trans_names, GENEID = GENEID)

# get all files
# files = Sys.glob( sprintf('%s/*/abundance.h5',input_dir) )
# files = sort( files )

# files_names = gsub('/abundance.h5','',files)
# files_names = gsub('.*/','',files_names)

files = sprintf('%s/%s',input_dir, df$SampleID)

df$files = files
df$files = sprintf('%s/abundance.h5',df$files)
df = df[order(df$files),]

# run tximport to summarize transcript to gene level
txi <- tximport(df$files, type = "kallisto", tx2gene = tx2gene)

temp = as.data.frame(txi$abundance); names(temp) = df$SampleID
temp = data.frame(Gene = row.names(temp),temp)
write.table(temp, file = sprintf('res/%s_txi_res.txt',name),quote = F, sep='\t',row.names = F)

# temp = as.data.frame(txi$counts); names(temp) = df$new_name
# temp = data.frame(Gene = row.names(temp),temp)
# write.table(temp, file = sprintf('res/%s_txi_counts_res_v2.txt',name),quote = F, sep='\t',row.names = F)

# filtering tximport results to wanted genes (protein coding from the table yael send, TPM1>0.2% samples)
gene_file = '/pita/users/tzipi/projects/rnaSeq/lncRNA/tables/from_yael/19815_Protein_coding_genes_GencodeV24.txt'
gene_df = read.table(gene_file, header = T)
wanted_genes = sprintf('%s_%s',gene_df$Gene,gene_df$Ensembl_Gene_ID) 
txi_genes = row.names(txi$abundance)
wanted_genes_pos1 = txi_genes %in% wanted_genes


## filter TPM to 1 TPM>0.2%
tpm_genes = row.names(txi$abundance)
wanted_genes_pos2 = vector(mode = 'logical',length = length(tpm_genes))
for ( i in 1:length(tpm_genes) )
{
  per = sum(txi$abundance[i,] > 1) / dim(txi$abundance)[2]
  wanted_genes_pos2[i] = per > 0.2
}
wanted_genes_pos = wanted_genes_pos1 & wanted_genes_pos2

txi$abundance = txi$abundance[wanted_genes_pos, ]
txi$counts = txi$counts[wanted_genes_pos, ]
txi$length = txi$length[wanted_genes_pos, ]
# 
name = sprintf('%s_geneFiltered', name)
temp = as.data.frame(txi$abundance); names(temp) = df$SampleID
temp = data.frame(Gene = row.names(temp),temp)
write.table(temp, file = sprintf('res/%s_txi_res.txt',name),quote = F, sep='\t',row.names = F)

# saveRDS(txi, file = sprintf('res/%s_txi_res.rds',name))
 

for (location_group in location_groups)
{
  txi_old = txi
  df_old = df
  name_old = name
  
  # filter to wanted locations and conditions
  name = sprintf('%s_%s_%s_vs_%s',name, location_group, two_groups[2],two_groups[1])
  pos = df[[prm]] %in% two_groups
  pos2 = df[[location_prm]] == location_group
  pos = pos & pos2
  
  if (one_sample_per_patient_filter_flag)
  {
    pos = pos[pos & df$one_sample_per_patient == 1]
    name = sprintf('%s_onePnSample', name)
  }

  df = df[pos,]
  txi$abundance = txi$abundance[, pos]
  txi$counts = txi$counts[, pos]
  txi$length = txi$length[, pos]
  

  coldata = data.frame(row.names=df$SampleID, diagnosis=as.factor(df[[prm]]) )
  
  ddsHTSeq <- DESeqDataSetFromTximport(txi = txi, colData = coldata, design = ~ diagnosis)
  colData(ddsHTSeq)$diagnosis<-factor(colData(ddsHTSeq)$diagnosis, levels=c(two_groups[1],two_groups[2]))
  
  
  dds<-DESeq(ddsHTSeq)
  
  res<-results(dds)
  res<-res[order(res$padj),]
  head(res)
  #
  # dir.create('res')
  out_path = 'DE/gene_filter_on_full_cohort/'
  write_DE_result(dds, name = name, out_path = out_path, type_flag = F)
  de_res = write_DE_result(dds, name = name, out_path = out_path, filter_flag = T, FDR_cutoff = 0.05, FC_cutoff = log2(1.5), type_flag = F)
  
  txi = txi_old
  df = df_old
  name = name_old
}

