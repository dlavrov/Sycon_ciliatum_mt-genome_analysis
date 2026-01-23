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
module load bowtie/1.3.0-jfngr2c
module load samtools/1.19.2-py311-qlmoxzv
module load bedtools2/2.31.1-py311-6kemgt3

# OUTPUT
output="$(basename -s .fasta ${input})-sam.out"

# STEP 1: Build Bowtie1 index
echo "Building Bowtie1 index..."
bowtie-build "$input" BOWTIELIB1

# STEP 2: Align reads strictly
echo "Running Bowtie1 alignment with strict settings..."
#bowtie -k 1 -p 4 --best BOWTIELIB1 --suppress 1,7 ${2},${3} bowtie1_out.sam
bowtie \
  -1 "${2}" -2 "${3}" \
  -a \
  -v 2 \ # allow only two mismatches, default -n 2 -l 28 two mismatches in the seed
  --best \ # reports best-scoring
  --strata \ #report only alignments with fewest mismatches 
  -X 1000 \ # this defines the fragement size for --fr option
  --fr \ # both forward and reverse should match (default option)
  -p "$threads" \
  -x BOWTIELIB1 \
  -S bowtie1_out.sam

# STEP 3: Convert to BAM, sort, and index
echo "Converting SAM to sorted BAM..."
samtools view -bS bowtie1_out.sam > bowtie1_out.bam
samtools sort bowtie1_out.bam -o bowtie1.sorted.bam
samtools index bowtie1.sorted.bam

# STEP 4: Create coverage file
echo "Calculating base-level coverage..."
samtools depth bowtie1.sorted.bam > bowtie1_coverage.txt

# STEP 5a: Average coverage per contig
echo "Calculating mean contig coverage..."
awk '{cov[$1]+=$3; len[$1]++} END {for (c in cov) print c, cov[c]/len[c]}' bowtie1_coverage.txt \
    | sort -k1.7n \
    > bowtie1_reads_mean_coverage.tsv

echo "✅ Alignment and coverage pipeline complete."
echo "➡ Output: bowtie1_reads_mean_coverage.tsv"

