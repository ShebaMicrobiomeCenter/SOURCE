# install.packages("rjson")
library("rjson")

metadata_file =  'data/rnaSeq_SOURCE_china_map.txt' 
metadata_df = read.table(file = metadata_file, header = T,sep = '\t', stringsAsFactors = F, na.strings = c('NA','na'))

metadata_df$fastq_read_count = vector(mode = 'numeric',length = dim(metadata_df)[1] )
for ( i in 1:dim(metadata_df)[1] )
{
  file_path = sprintf('kallisto/%s/run_info.json', metadata_df$SampleID[i])
  json_data <- fromJSON(paste(readLines(file_path), collapse=""))
  metadata_df$fastq_read_count[i] = json_data$n_processed
}

metadata_df$post_kallisto_read_count = vector(mode = 'numeric',length = dim(metadata_df)[1] )
for ( i in 1:dim(metadata_df)[1] )
{
  file_path = sprintf('kallisto/%s/abundance.tsv', metadata_df$SampleID[i])
  df = read.table(file = file_path, header = T,sep = '\t')
  metadata_df$post_kallisto_read_count[i] = sum(df$est_counts)
  # print(i)
}


metadata_df = metadata_df[, c(1,( dim(metadata_df)[2]-1 ), ( dim(metadata_df)[2] ),  2:( dim(metadata_df)[2]-2 ) ) ]

out_file = 'data/rnaSeq_SOURCE_china_map_reads.txt' 
write.table(x = metadata_df, file = out_file, quote = F, sep='\t', row.names = F)
