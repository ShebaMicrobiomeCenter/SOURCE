## Running KneadData version v0.12.0:
# Running KneadData for paired end samples (Chinese SOURCE):
for f in /PATH/TO/RAW/DATA/*_R1.fastq; do f=${f##*/}; f=${f%%_R*};  kneaddata --input1 /PATH/TO/RAW/DATA/${f}_R1.fastq --input2 /PATH/TO/RAW/DATA/${f}_R2.fastq --output /OUTPUT/DIR/kneaddata/main --scratch /OUTPUT/DIR/tmp --threads 70 --output-prefix ${f}  --cat-final-output --reference-db /REFERENCE/DIR/kneaddata_db/human_genome_bowtie2 --trimmomatic /PATH/TO/CONDA/miniconda3/envs/biobakery3.1/share/trimmomatic-0.39-2/  --serial --run-trf --remove-intermediate-output; done

# Running KneadData for single end samples (Israeli SOURCE):
for f in /PATH/TO/RAW/DATA2/*_R1.fastq.gz
do
f=${f##*/}
f=${f%%_R*}
kneaddata --un /PATH/TO/RAW/DATA2/${f}_R1.fastq.gz --output /OUTPUT/DIR2/kneaddata/main --threads 70 --output-prefix ${f} --reference-db /REFERENCE/DIR/human_genome_bowtie2 --trimmomatic /PATH/TO/CONDA/miniconda3/envs/huttenhower-py3.7/share/trimmomatic-0.39-2/  --serial --run-trf  --remove-intermediate-output
done

# Combining KneadData results:
kneaddata_read_count_table --input /OUTPUT/DIR/kneaddata/main --output /OUTPUT/DIR/kneaddata/merged/kneaddata_read_count_table.tsv

## Running MetaPhlAn version v4.0.0 on Chinese SOURCE:
for f in /PATH/TP/RAW/DATA/*_R1.fastq
do
f=${f##*/}
f=${f%%_R*}
metaphlan /OUTPUT/DIR/kneaddata/main/${f}.fastq --input_type fastq --output_file /OUTPUT/DIR/metaphlan/main/${f}_taxonomic_profile.tsv --samout /OUTPUT/DIR/metaphlan/main/${f}_bowtie2.sam --nproc 50 --no_map --tmp_dir /OUTPUT/DIR/metaphlan/main --bowtie2db=/REFERENCE/DIR/metaphlan_db/v31 --index=mpa_vJan21_CHOCOPhlAnSGB_202103
done

# Merging MetaPhlAn abundance tables of all samples into one table: (this is the file used for analysis)
cd /OUTPUT/DIR/metaphlan/main:
merge_metaphlan_tables.py *_taxonomic_profile.tsv > ../merged/metaphlan_taxonomic_profiles.tsv

## Running HUMAnN version v3.6 on Chinese SOURCE:
for f in /PATH/TO/RAW/DATA/*_R1.fastq
do
f=${f##*/}
f=${f%%_R*}
humann -i /OUTPUT/DIR/kneaddata/main/${f}.fastq -o /OUTPUT/DIR/humann/main --o-log /OUTPUT/DIR/humann/main/${f}.log --threads 70 --taxonomic-profile /OUTPUT/DIR/metaphlan/main/${f}_taxonomic_profile.tsv --input-format fastq --remove-temp-output
done

# Joining HUMAnN tables for different features: pathabundance, and genefamilies (switch 'pathabundance' with 'genefamilies'):
humann_join_tables --input /OUTPUT/DIR/humann/main/ --output /OUTPUT/DIR/humann/merged/pathabundance.tsv --file_name pathabundance

# ecs tables requires regrouping from the genefamilies abundance table:
humann_regroup_table --input /OUTPUT/DIR/humann/merged/genefamilies.tsv --output /OUTPUT/DIR/humann/merged/ecs.tsv --groups uniref90_level4ec

# Renorming pathabundance/genefamilies/ecs tables:
humann_renorm_table --input /OUTPUT/DIR/humann/merged/pathabundance.tsv --output /OUTPUT/DIR/humann/relab/pathabundance_relab.tsv --units relab --special n

# Splitting feature abundance tables into stratified tables (i.e. tables of features and their origin species) and unstratified tables:
humann_split_stratified_table -i /OUTPUT/DIR/relab/pathabundance_relab.tsv -o /OUTPUT/DIR/humann/relab/

# Files used for analysis are the unstratified abundance tables which were outputted by the command humann_split_stratified_table, e.g., pathabundance_relab_unstratified.tsv
