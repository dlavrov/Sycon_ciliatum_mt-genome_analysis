#!/bin/bash

# Usage check
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 input_multifasta.fa output_results.txt"
  exit 1
fi

input_fasta=$1
output_file=$2
temp_fasta="temp_single_seq.fa"

# Clear output file if exists
> "$output_file"

# Extract sequence IDs
#seq_ids=$(grep "^>" "$input_fasta" | sed 's/^>//')

echo "Processing sequences from $input_fasta..."

grep "^>" "$input_fasta" | while read -r header; do
  id=${header#>}       # remove leading >
  id=${id%% *}         # keep only first token

  echo "Processing sequence: $id"

  awk -v seqid=">$id" '
    BEGIN {flag=0}
    $0 ~ "^" seqid {flag=1; print; next}
    /^>/ && flag {exit}
    flag {print}
  ' "$input_fasta" > "$temp_fasta"

  {
    echo "### Results for sequence: $id ###"
    etandem -sequence "$temp_fasta" -outfile stdout -minrepeat 20 -maxrepeat 200
    echo
  } >> "$output_file"

done

# Cleanup
rm -f "$temp_fasta"

echo "All done. Results compiled in $output_file"
