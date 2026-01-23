#!/usr/bin/env python3

import sys
import re
from Bio import SeqIO

def parse_note_field(note):
    """Parse Note= field in GFF to extract size and copy number."""
    size_match = re.search(r"Size:(\d+)", note)
    copies_match = re.search(r"Copies:(\d+)", note)
    size = int(size_match.group(1)) if size_match else None
    copies = int(copies_match.group(1)) if copies_match else None
    return size, copies

def extract_repeats(gff_file, fasta_file, min_size, max_size):
    """Extract tandem repeats within the given size range and save to a FASTA file."""
    # Load FASTA
    fasta_seqs = {record.id: record.seq for record in SeqIO.parse(fasta_file, "fasta")}
   
    # Construct output filename automatically
    output_fasta = f"etandem_{min_size}-{max_size}_repeats.fa"

    # Prepare output
    with open(output_fasta, "w") as out_f:
        with open(gff_file) as gff:
            for line in gff:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.strip().split("\t")
                if len(fields) != 9:
                    continue  # skip malformed lines

                chrom = fields[0].strip().split()[0]  # remove anything after space (e.g., '###')
                start, end = int(fields[3]), int(fields[4])
                attributes = fields[8]

                note = attributes.split("Note=")[1] if "Note=" in attributes else ""
                size, copies = parse_note_field(note)
                if size is None or copies is None:
                    continue

                print(f"Checking {chrom}: size={size}, count={copies}, start={start}, end={end}")
                if not (min_size <= size <= max_size):
                    continue

                sequence = fasta_seqs.get(chrom)
                if not sequence:
                    print(f"WARNING: Chromosome {chrom} not found in FASTA")
                    continue

                repeat_region = sequence[start - 1:end]  # GFF is 1-based
                expected_len = size * copies
                if len(repeat_region) < expected_len:
                    print(f"WARNING: Region too short for {copies} copies of size {size}")
                    continue

                # Split into individual repeats
                for i in range(copies):
                    print(f"Writing chunk {i + 1} of {copies} for {chrom}_{start}_{end}")
                    repeat = repeat_region[i * size : (i + 1) * size]
                    print(f"Expected size: {size}, actual repeat length: {len(repeat)}")
                    if len(repeat) != size:
                        continue  # skip partial repeat
                    header = f">{chrom}_{start}_{end}_{i + 1}"
                    out_f.write(f"{header}\n{repeat}\n")

    print(f"Repeats written to {output_fasta}")

# Example usage:
# python extract_repeats.py mtseq.1.gff mtseq.1.fasta repeats_out.fasta

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python extract_repeats.py input.gff input.fasta output.fasta [min_size] [max_size")
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

