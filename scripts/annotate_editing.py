#!/usr/bin/env python3
import sys
import argparse
from textwrap import wrap

# ---------------------------------------------------------
# CLI
# ---------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Annotate and edit mitochondrial RNA editing sites in FASTA"
)
parser.add_argument("input", help="Input FASTA")
parser.add_argument("output", help="Output FASTA")
parser.add_argument("--mode", choices=["inline", "gff"], default="inline",
                    help="Annotation mode (default: inline)")
parser.add_argument("--gff", help="GFF3 output file (required if --mode gff)")
parser.add_argument("--wrap", type=int, default=0,
                    help="FASTA line width (0 = no wrapping)")
args = parser.parse_args()

if args.mode == "gff" and not args.gff:
    sys.exit("ERROR: --gff is required when --mode gff is used")

# ---------------------------------------------------------
# Annotation helpers
# ---------------------------------------------------------
def inline_annotation(type_no, idx, motif):
    return [
        f";     editing_site_type{type_no}_{idx} ==> start",
        motif,
        f";     editing_site_type{type_no}_{idx} ==> end"
    ]

def gff_entry(seqid, start, end, type_no, idx, inserted, context):
    attrs = (
        f"ID=editing_site_type{type_no}_{idx};"
        f"Type=type{type_no};"
        f"InsertedTs={inserted};"
        f"Context={context}"
    )
    return (
        f"{seqid}\tmt_editing\tediting_site_type{type_no}\t"
        f"{start}\t{end}\t.\t+\t.\t{attrs}"
    )

# ---------------------------------------------------------
# Core processing
# ---------------------------------------------------------
def process_sequence(seq, seqid, wrap_width, mode, gff_out):
    edited_seq = []
    inline_lines = []
    buf = []
    counters = {1: 0, 2: 0, 3: 0}

    def flush_buf():
        nonlocal buf
        if buf:
            edited_seq.append("".join(buf))
            buf = []

    i = 0
    while i < len(seq):

        if seq[i] == "T":
            j = i
            while j < len(seq) and seq[j] == "T":
                j += 1

            run_len = j - i
            type_no = None
            insert = 0

            # ≥4 Ts → type 3
            if run_len >= 4:
                type_no, insert = 3, 2

            # exactly 3 Ts
            elif run_len == 3 and i > 0:
                prev = seq[i - 1]
                if prev in "ACT":
                    type_no, insert = 3, 2
                elif prev == "G":
                    type_no, insert = 2, 1

            # exactly 2 Ts with CCC upstream
            elif run_len == 2 and i >= 3 and seq[i - 3:i] == "CCC":
                type_no, insert = 1, 1

            if type_no:
                flush_buf()

                counters[type_no] += 1
                edited = "T" * run_len + "t" * insert
                edited_seq.append(edited)

                start = i + 1
                end = i + run_len + insert

                if mode == "inline":
                    inline_lines.extend(
                        inline_annotation(type_no, counters[type_no], edited)
                    )
                else:
                    context = seq[max(0, i - 3):j]
                    gff_out.append(
                        gff_entry(
                            seqid, start, end,
                            type_no, counters[type_no],
                            insert, context
                        )
                    )

                i = j
                continue

            # not an editing site
            buf.append(seq[i])
            i += 1
            continue

        buf.append(seq[i])
        i += 1

    flush_buf()

    full_edited = "".join(edited_seq)
    fasta_lines = wrap(full_edited, wrap_width) if wrap_width else [full_edited]

    return fasta_lines, inline_lines

# ---------------------------------------------------------
# FASTA I/O
# ---------------------------------------------------------
gff_lines = ["##gff-version 3"]

with open(args.input) as fin, open(args.output, "w") as fout:
    header = None
    seq = []

    for line in fin:
        line = line.rstrip()
        if line.startswith(">"):
            if header:
                seqid = header[1:].split()[0]
                fasta_lines, inline_lines = process_sequence(
                    "".join(seq), seqid, args.wrap, args.mode, gff_lines
                )

                fout.write(header + "\n")
                for l in fasta_lines:
                    fout.write(l + "\n")

                if args.mode == "inline":
                    for l in inline_lines:
                        fout.write(l + "\n")

            header = line
            seq = []
        else:
            seq.append(line)

    # last record
    if header:
        seqid = header[1:].split()[0]
        fasta_lines, inline_lines = process_sequence(
            "".join(seq), seqid, args.wrap, args.mode, gff_lines
        )

        fout.write(header + "\n")
        for l in fasta_lines:
            fout.write(l + "\n")

        if args.mode == "inline":
            for l in inline_lines:
                fout.write(l + "\n")

# ---------------------------------------------------------
# Write GFF if requested
# ---------------------------------------------------------
if args.mode == "gff":
    with open(args.gff, "w") as gff:
        gff.write("\n".join(gff_lines) + "\n")

print("Done.")

