#!/bin/bash
set -ex #exit as soon as any line fails and print each command

if [ "$#" -lt 1 ] 
then 
	echo "error: a name of a SRA file is required"
	exit 1
fi

############################ SPLIT SRA ##############################
 
module load sratoolkit
fastq-dump --split-files $1

# #Make some variables
# 
# #filename=`echo "$1" | cut -d'.' -f1`
filename="${1%.*}"
# 
# #can be done like this also:
# #filename=$(basename "$fullfile")
# #extension="${filename##*.}"
# #filename="${filename%.*}"
# 
#this can be replaced with files=`ls *.fastq`
# 
[ -d raw ] && echo "Directory 'raw' exists" || mkdir raw
mv $filename'_1.fastq' raw/$filename'_1.fq'
mv $filename'_2.fastq' raw/$filename'_2.fq'

############################ FASTQC #######################
# 
[ -d qc ] && echo "Directory 'qc' exists" || mkdir qc
[ -d qc/raw ] && echo "Directory 'qc/raw' exists" || mkdir qc/raw
# 
 module load fastqc/
 fastqc -o qc/raw/ -t 8 raw/*.fq
# 
########################### SortMeRNA ###########################
# 
# smr_home='/home/dlavrov/bin/sortmerna-2.1-linux-64'
# export SORTMERNA_DB="$smr_home/rRNA_databases/silva-bac-16s-id90.fasta,$smr_home/index/silva-bac-16s-db:$smr_home/rRNA_databases/silva-bac-23s-id98.fasta,$smr_home/index/silva-bac-23s-db:$smr_home/rRNA_databases/silva-arc-16s-id95.fasta,$smr_home/index/silva-arc-16s-db:$smr_home/rRNA_databases/silva-arc-23s-id98.fasta,$smr_home/index/silva-arc-23s-db:$smr_home/rRNA_databases/silva-euk-18s-id95.fasta,$smr_home/index/silva-euk-18s-db:$smr_home/rRNA_databases/silva-euk-28s-id98.fasta,$smr_home/index/silva-euk-28s:$smr_home/rRNA_databases/rfam-5s-database-id98.fasta,$smr_home/index/rfam-5s-db:$smr_home/rRNA_databases/rfam-5.8s-database-id98.fasta,$smr_home/index/rfam-5.8s-db"
#module load sortmerna/
# $smr_home/scripts/merge-paired-reads.sh raw/$filename'_1.fq' raw/$filename'_2.fq' $filename'_merged.fq'
# $smr_home/sortmerna --ref $SORTMERNA_DB --reads $filename'_merged.fq' --paired_in -a 16 --log --fastx --aligned $filename'_rRNA' --other $filename'_sortmerna'
# $smr_home/scripts/unmerge-paired-reads.sh $filename'_sortmerna.fq' $filename'_sortmerna_1.fq' $filename'_sortmerna_2.fq'
# rm *merged.fq
# rm *_sortmerna.fq
# mkdir sortmerna
# mv *sortmerna*fq sortmerna
# mv *rRNA* sortmerna
# filename=$filename'_sortmerna'
# 
############################ FASTQC #################################
# 
# mkdir qc/sortmerna
# fastqc -o qc/sortmerna/ -t 8 sortmerna/$filename*.fq
# 
# ############################ TRIMMOMATIC ############################

db_home='/home/dlavrov/share/trimmomatic'
module unload java/1.7.0_55
module load trimmomatic
trimmomatic  PE -threads 8 raw/$filename'_1.fq' raw/$filename'_2.fq' $filename'_trimmomatic_1.fq' $filename'_trimmomatic_u1.fq' $filename'_trimmomatic_2.fq' $filename'_trimmomatic_u2.fq'  ILLUMINACLIP:$db_home/TruSeq3-PE-2.fa:3:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50

[ -d trimmomatic ] && echo "Directory 'trimmomatic' xists" || mkdir trimmomatic
mv *_trimmomatic_* trimmomatic
[ -d qc/trimmomatic ] && echo "Directory 'qc/trimmomatic' exists" || mkdir qc/trimmomatic
filename=$filename'_trimmomatic'

############################ FASTQC #################################

module load fastqc
fastqc -o qc/trimmomatic -t 8 trimmomatic/$filename*.fq

###################### Renaming sequences ###########################
# 
# sed -i -r 's/\s+/_/g' trimmomatic/$filename'_1.fq'
# sed -i -r 's/_length/\\1 length/' trimmomatic/$filename'_1.fq'
# sed -i -r 's/\s+/_/g' trimmomatic/$filename'_2.fq'
# sed -i -r 's/_length/\\2 length/' trimmomatic/$filename'_2.fq'

