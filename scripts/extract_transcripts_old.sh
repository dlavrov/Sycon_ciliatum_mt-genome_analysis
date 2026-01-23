#!/bin/bash
# Usage: ./extract_transcripts.sh blast_output.txt gene_map.csv

# --- Argument check ---
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <blast_output.txt> <gene_map.csv>"
    echo
    echo "  blast_output.txt : BLAST output file with at least 4 columns"
    echo "  gene_map.csv     : CSV mapping file with two columns (e.g., KU244283.1,atp6)"
    exit 1
fi

blast_file=$1
map_file=$2

# --- Check input files exist ---
if [ ! -f "$blast_file" ]; then
    echo "Error: BLAST file '$blast_file' not found."
    exit 1
fi

if [ ! -f "$map_file" ]; then
    echo "Error: mapping file '$map_file' not found."
    exit 1
fi

# --- Main processing ---
awk -F'[ ,\t]+' '
    BEGIN { OFS="\t" }   # Output fields separated by tabs

    # Phase 1: read mapping file (CSV)
    FNR==NR {
        id2gene[$1] = $2
        next
    }

    # Phase 2: process BLAST file
    {
        id = $1
        gene = id2gene[id]
        if (gene != "") {
            outfile = gene "_transcript.lst"
            print $2, $3, $4 >> outfile
        }
    }
' "$map_file" "$blast_file"

echo "Processing complete."

