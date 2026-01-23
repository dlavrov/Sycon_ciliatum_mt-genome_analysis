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
  -1 "${reads_R1}" -2 "${reads_R2}" \
  -a \
  -v 2 \
  --best \
  --strata \
  --suppress 1,7 \
  -X 1000 \
  --fr \
  -p "$threads" \
  -x BOWTIELIB1 \
  bowtie1_out.tab

