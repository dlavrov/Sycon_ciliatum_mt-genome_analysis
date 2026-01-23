#!/usr/bin/env python3
import argparse
import sys
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(
        description="Extract feature sequences from genome based on BLAST-like coordinates."
    )
    parser.add_argument("genome", help="Genome FASTA file")
    parser.add_argument("hits", help="BLAST results file with columns: seq_name start end")
    parser.add_argument(
        "-o", "--output", help="Output FASTA file (default: print to stdout)"
    )
    args = parser.parse_args()

    # Load genome sequences into dictionary
    genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))

    # Determine output destination
    out_handle = open(args.output, "w") if args.output else sys.stdout

    with open(args.hits) as hits, out_handle:
        for line in hits:
            if not line.strip() or line.startswith("#"):
                continue
            name, start_str, end_str = line.strip().split("\t")
            start, end = int(start_str), int(end_str)

            if name not in genome:
                raise ValueError(f"Sequence {name} not found in genome FASTA")

            seq_record = genome[name]

            if start < end:
                # forward orientation
                subseq = seq_record.seq[start - 1:end]
                strand = "f"
            else:
                # reverse orientation
                subseq = seq_record.seq[end - 1:start].reverse_complement()
                strand = "r"

            header = f"{name}_{min(start, end)}_{max(start, end)}_{strand}"
            print(f">{header}\n{subseq}", file=out_handle)

    # Only print message if writing to file
    if args.output:
        print(f"Extracted features saved to {args.output}")

if __name__ == "__main__":
    main()

