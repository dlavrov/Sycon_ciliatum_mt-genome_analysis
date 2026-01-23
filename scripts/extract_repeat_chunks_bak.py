#!/usr/bin/env python3

import sys
import re
from Bio import SeqIO

def parse_gff_attributes(attr_string):
    """Parse GFF3 attributes into a dictionary."""
    attrs = {}
    for field in attr_string.split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            attrs[key] = value
    return attrs

def extract_repeats(gff_file, fasta_file, min_size, max_size):
    """Extract tandem repeats within the given size range and save to a FASTA file."""

    # Load FASTA
    # fasta_seqs = {rec.id: rec.seq for rec in SeqIO.parse(fasta_file, "fasta")}
    fasta_seqs = {
        record.description.split()[0]: record.seq
        for record in SeqIO.parse(fasta_file, "fasta")
    }

    output_fasta = f"etandem_{min_size}-{max_size}_repeats.fa"

    with open(output_fasta, "w") as out_f, open(gff_file) as gff:
        for line in gff:
            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip().split("\t")
            if len(fields) != 9:
                continue

            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            attributes = parse_gff_attributes(fields[8])

            # Required attributes
            if "Size" not in attributes or "Copies" not in attributes:
                continue

            size = int(attributes["Size"])
            copies = int(attributes["Copies"])
            feature_id = attributes.get("ID", f"{chrom}_{start}_{end}")

            if not (min_size <= size <= max_size):
                continue

            sequence = fasta_seqs.get(chrom)
            if sequence is None:
                print(f"WARNING: {chrom} not found in FASTA")
                continue

            repeat_region = sequence[start - 1:end]  # GFF is 1-based
            expected_len = size * copies

            if len(repeat_region) < expected_len:
                print(f"WARNING: {feature_id} region too short")
                continue

            # Split into individual repeat units
            for i in range(copies):
                repeat = repeat_region[i * size:(i + 1) * size]
                if len(repeat) != size:
                    continue

                header = f">{feature_id}_copy{i + 1}"
                out_f.write(f"{header}\n{repeat}\n")

    print(f"Repeats written to {output_fasta}")

# ------------------------------------------------------------
# CLI
# ------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: extract_repeats.py input.gff input.fasta [min_size] [max_size]")
        sys.exit(1)

    gff_file = sys.argv[1]
    fasta_file = sys.argv[2]

    if len(sys.argv) == 5:
        min_size = int(sys.argv[3])
        max_size = int(sys.argv[4])
    else:
        min_size = 10
        max_size = 200

    extract_repeats(gff_file, fasta_file, min_size, max_size)

