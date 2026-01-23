#!/usr/bin/env python3
import argparse
from Bio import SeqIO

def find_runs(seq, base, min_len):
    """Return list of (start1, end1) 1-based inclusive runs of base in seq (case-insensitive)."""
    seqU = str(seq).upper()
    runs = []
    i = 0
    n = len(seqU)
    while i < n:
        if seqU[i] == base:
            j = i
            while j < n and seqU[j] == base:
                j += 1
            run_len = j - i
            if run_len >= min_len:
                runs.append((i + 1, j))  # 1-based inclusive
            i = j
        else:
            i += 1
    return runs

def main(fasta, out_gff, min_len, end_window):
    with open(out_gff, "w") as out:
        out.write("##gff-version 3\n")
        out.write("# polyA/polyT mononucleotide runs (A-only or T-only)\n")

        for rec in SeqIO.parse(fasta, "fasta"):
            seqid = rec.id
            L = len(rec.seq)

            for base, ftype in [("A", "polyA_run"), ("T", "polyT_run")]:
                runs = find_runs(rec.seq, base, min_len)
                idx = 0
                for (s, e) in runs:
                    # If end_window is set, keep only runs within end_window of either end
                    if end_window is not None:
                        near_left = s <= end_window
                        near_right = e >= (L - end_window + 1)
                        if not (near_left or near_right):
                            continue

                    idx += 1
                    # Mark which end it is closest to (helpful for plotting)
                    if end_window is not None:
                        end_tag = "left" if s <= end_window else ("right" if e >= (L - end_window + 1) else "internal")
                    else:
                        # still useful:
                        end_tag = "left" if s <= min_len else ("right" if e >= L - min_len + 1 else "internal")

                    attrs = (
                        f"ID={ftype}_{seqid}_{idx};"
                        f"Base={base};Length={e-s+1};End={end_tag}"
                    )
                    out.write(f"{seqid}\tpolyAT\t{ftype}\t{s}\t{e}\t.\t+\t.\t{attrs}\n")

    print(f"Wrote {out_gff}")

if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Detect mononucleotide A-only / T-only runs in FASTA and output GFF3."
    )
    p.add_argument("fasta", help="Input FASTA (contigs/chromosomes)")
    p.add_argument("gff", help="Output GFF3")
    p.add_argument("--min-len", type=int, default=20,
                   help="Minimum run length to report (default: 20)")
    p.add_argument("--end-window", type=int, default=200,
                   help="Only report runs within N bp of either contig end (default: 200). "
                        "Use 0 or omit to report anywhere.")
    args = p.parse_args()

    ew = args.end_window
    if ew == 0:
        ew = None

    main(args.fasta, args.gff, args.min_len, ew)

