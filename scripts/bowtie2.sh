#!/bin/bash
set -euo pipefail

#############################################
#           run_alignment_pipeline.sh       #
#############################################

# USAGE CHECK
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 contigs_file short_reads_file1 short_reads_file2"
    exit 1
fi

# INPUTS
threads=32
input=$1      # contigs file (FASTA)
reads_R1=$2   # paired-end reads 1
reads_R2=$3   # paired-end reads 2

# MODULES
module load samtools/1.19.2-py311-qlmoxzv
module load bedtools2/2.31.1-py311-6kemgt3

# STEP 1: Build Bowtie2 index
echo "Building Bowtie2 index..."
bowtie2-build "$input" BOWTIELIB2

# STEP 2: Align reads strictly
echo "Running Bowtie2 alignment with strict settings..."
bowtie2 -x BOWTIELIB2 \
        -1 "$reads_R1" -2 "$reads_R2" \
        --end-to-end \
        --no-unal \
        --no-mixed \
        --no-discordant \
        -k 1 \
        -p "$threads" \
        -S aligned.sam

# STEP 3: Convert to BAM, sort, and index
echo "Converting SAM to sorted BAM..."
samtools view -bS aligned.sam > aligned.bam
samtools sort aligned.bam -o aligned.sorted.bam
samtools index aligned.sorted.bam

# STEP 4: Create coverage file
echo "Calculating base-level coverage..."
samtools depth aligned.sorted.bam > coverage.txt

# STEP 5a: Average coverage per contig
echo "Calculating mean contig coverage..."
awk '{cov[$1]+=$3; len[$1]++} END {for (c in cov) print c, cov[c]/len[c]}' coverage.txt \
    | sort -k1.7n \
    > reads_mean_coverage.tsv

echo "✅ Alignment and coverage pipeline complete."
echo "➡ Output: reads_mean_coverage.tsv"

