#!/usr/bin/env python3

import sys
import re
import argparse
from Bio import SeqIO


def parse_note_field(note):
    """Parse Note= field in GFF to extract size and copy number."""
    size_match = re.search(r"Size[:=](\d+)", note)
    copies_match = re.search(r"Copies[:=](\d+)", note)
    size = int(size_match.group(1)) if size_match else None
    copies = int(copies_match.group(1)) if copies_match else None
    return size, copies


def extract_repeats(gff_file, fasta_file, min_size, max_size, verbose=False):
    """Extract tandem repeats within the given size range and save to a FASTA file."""

    # Load FASTA
    fasta_seqs = {record.id: record.seq for record in SeqIO.parse(fasta_file, "fasta")}

    output_fasta = f"etandem_{min_size}-{max_size}_repeats.fa"

    with open(output_fasta, "w") as out_f, open(gff_file) as gff:
        for line in gff:
            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip().split("\t")
            if len(fields) != 9:
                continue

            chrom = fields[0].split()[0]
            start, end = int(fields[3]), int(fields[4])
            attrs = fields[8]

            size, copies = parse_note_field(attrs)
            if size is None or copies is None:
                continue

            if not (min_size <= size <= max_size):
                continue

            sequence = fasta_seqs.get(chrom)
            if sequence is None:
                print(f"WARNING: Chromosome {chrom} not found in FASTA", file=sys.stderr)
                continue

            repeat_region = sequence[start - 1:end]
            expected_len = size * copies

            if len(repeat_region) < expected_len:
                print(
                    f"WARNING: Region too short for {copies} copies of size {size} "
                    f"({chrom}:{start}-{end})",
                    file=sys.stderr
                )
                continue

            for i in range(copies):
                repeat = repeat_region[i * size:(i + 1) * size]

                if verbose:
                    print(
                        f"Writing chunk {i + 1} of {copies} "
                        f"for {chrom}_{start}_{end}"
                    )
                    print(
                        f"Expected size: {size}, actual repeat length: {len(repeat)}"
                    )

                if len(repeat) != size:
                    continue

                header = f">{chrom}_{start}_{end}_{i + 1}"
                out_f.write(f"{header}\n{repeat}\n")

    if verbose:
        print(f"Repeats written to {output_fasta}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract individual tandem repeat units from GFF + FASTA"
    )
    parser.add_argument("gff", help="Input GFF file with tandem repeats")
    parser.add_argument("fasta", help="Genome FASTA file")
    parser.add_argument("min_size", nargs="?", type=int, default=10,
                        help="Minimum repeat unit size (default: 10)")
    parser.add_argument("max_size", nargs="?", type=int, default=200,
                        help="Maximum repeat unit size (default: 200)")
    parser.add_argument("--verbose", action="store_true",
                        help="Print progress and diagnostic information")

    args = parser.parse_args()

    extract_repeats(
        args.gff,
        args.fasta,
        args.min_size,
        args.max_size,
        verbose=args.verbose
    )


if __name__ == "__main__":
    main()

