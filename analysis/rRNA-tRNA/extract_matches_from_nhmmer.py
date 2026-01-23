#!/usr/bin/env python3

import argparse
import sys
from Bio import SeqIO
from Bio.Seq import Seq


def wrap_fasta(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main():
    p = argparse.ArgumentParser(
        description=(
            "Extract regions from FASTA using a whitespace-delimited table.\n"
            "Uses: seqid=field1, start=field7, end=field8, strand=field12.\n"
            "Reverse-complements extracted sequence if strand is '-'."
        )
    )
    p.add_argument("table", help="Input table file (like nhmmer tblout-style; comments start with #)")
    p.add_argument("fasta", help="FASTA file with target sequences")
    p.add_argument("-o", "--out", default="extracted.fa", help="Output FASTA (default: extracted.fa)")
    p.add_argument("--minlen", type=int, default=1, help="Skip extracted regions shorter than this (default: 1)")
    p.add_argument("--keep-original-order", action="store_true",
                   help="Do not normalize coordinates; use start..end as given (rarely desired).")
    args = p.parse_args()

    # Load FASTA into memory (fast for typical mt contig sets; for huge genomes you’d stream)
    seqs = {rec.id: rec.seq for rec in SeqIO.parse(args.fasta, "fasta")}

    n_written = 0
    n_skipped = 0

    with open(args.table) as fin, open(args.out, "w") as fout:
        for line_no, line in enumerate(fin, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            fields = line.split()
            if len(fields) < 12:
                print(f"WARNING: line {line_no}: expected >=12 fields, got {len(fields)}; skipping", file=sys.stderr)
                n_skipped += 1
                continue

            seqid = fields[0]

            try:
                start_raw = int(fields[6])
                end_raw = int(fields[7])
            except ValueError:
                print(f"WARNING: line {line_no}: start/end not integers; skipping", file=sys.stderr)
                n_skipped += 1
                continue

            strand = fields[11]
            if strand not in {"+", "-"}:
                print(f"WARNING: line {line_no}: strand '{strand}' not '+' or '-'; skipping", file=sys.stderr)
                n_skipped += 1
                continue

            if seqid not in seqs:
                print(f"WARNING: line {line_no}: seq '{seqid}' not found in FASTA; skipping", file=sys.stderr)
                n_skipped += 1
                continue

            seq = seqs[seqid]

            if args.keep_original_order:
                # Use coords as given; still enforce bounds after converting to 0-based slices.
                start = start_raw
                end = end_raw
            else:
                # Normalize to start <= end for slicing.
                start = min(start_raw, end_raw)
                end = max(start_raw, end_raw)

            # Bounds check (1-based inclusive)
            if start < 1 or end > len(seq):
                print(
                    f"WARNING: line {line_no}: interval {seqid}:{start}-{end} out of bounds (len={len(seq)}); skipping",
                    file=sys.stderr
                )
                n_skipped += 1
                continue

            # Extract (convert to python slice: 0-based, end-exclusive)
            subseq = seq[start - 1:end]

            # Reverse-complement if on minus strand
            if strand == "-":
                subseq = subseq.reverse_complement()

            if len(subseq) < args.minlen:
                n_skipped += 1
                continue

            # Header includes normalized coords + strand; also record original raw coords for traceability
            header = f"{seqid}:{start_raw}-{end_raw}|norm:{start}-{end}|strand:{strand}"
            fout.write(f">{header}\n{wrap_fasta(str(subseq))}\n")
            n_written += 1

    print(f"Wrote {n_written} sequences to {args.out} (skipped {n_skipped}).", file=sys.stderr)


if __name__ == "__main__":
    main()

