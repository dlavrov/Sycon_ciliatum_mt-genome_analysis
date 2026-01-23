#!/usr/bin/env bash
# Extract transcript information from BLAST output and write to gene-specific files.
# Usage: ./extract_transcripts.sh blast_output.txt

# --- Check that one argument is provided ---
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 blast_output.txt"
    exit 1
fi

blast_file="$1"

# --- Define the mapping directly in the script ---
declare -A gene_map=(
    [KU244283.1]="atp6"
    [KU244282.1]="nad4"
    [KU244281.1]="nad5"
    [KU244280.1]="cox2"
    [KU244279.1]="cox3"
    [KU244278.1]="nad1"
    [KU244277.1]="nad3"
    [KU244276.1]="cox1"
    [KU244275.1]="cob"
    [KU244274.1]="nad2"
    [KU244273.1]="atp9"
)

# --- Remove any old *_transcript.lst files for these genes ---
for gene in "${gene_map[@]}"; do
    rm -f "${gene}_transcript.lst"
done

# --- Initialize a counter array ---
declare -A counts

# --- Process each line of the BLAST output ---
while read -r line; do
    # Skip empty lines or comments
    [[ -z "$line" || "$line" =~ ^# ]] && continue

    # Extract fields (allow tab or space delimiters)
    seq_id=$(echo "$line" | awk '{print $1}')
    field2=$(echo "$line" | awk '{print $2}')
    field3=$(echo "$line" | awk '{print $3}')
    field4=$(echo "$line" | awk '{print $4}')

    # Look up gene name
    gene_name="${gene_map[$seq_id]}"
    if [[ -z "$gene_name" ]]; then
        echo "Warning: No gene name found for $seq_id" >&2
        continue
    fi

    # Write tab-delimited line to gene-specific file
    printf "%s\t%s\t%s\n" "$field2" "$field3" "$field4" >> "${gene_name}_transcript.lst"

    # Increment counter
    ((counts["$gene_name"]++))
done < "$blast_file"

# --- Print summary ---
echo
echo "Transcript extraction completed."
echo "---------------------------------"
printf "%-10s %s\n" "Gene" "Entries"
echo "---------------------------------"

for gene in "${!gene_map[@]}"; do
    name="${gene_map[$gene]}"
    count="${counts[$name]:-0}"
    printf "%-10s %d\n" "$name" "$count"
done | sort

echo "---------------------------------"

