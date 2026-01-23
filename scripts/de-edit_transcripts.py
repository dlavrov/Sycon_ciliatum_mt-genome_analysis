#!/usr/bin/env python3
import argparse
from textwrap import wrap

# ---------------------------------------------------------
# CLI
# ---------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Reverse mitochondrial RNA editing (remove inserted Ts)"
)
parser.add_argument("input", help="Edited FASTA")
parser.add_argument("output", help="Unedited FASTA")
parser.add_argument("--wrap", type=int, default=0,
                    help="FASTA line width (0 = no wrapping)")
args = parser.parse_args()

# ---------------------------------------------------------
# Core unediting logic
# ---------------------------------------------------------
def unedit_sequence(seq):
    """
    Case-insensitive unediting rules:

    Type 3:
      - last 5 Ts of a run are preceded by A/C/T
      - remove 2 Ts from the end

    Type 2:
      - exactly 4 Ts
      - preceded by G
      - remove 1 T

    Type 1:
      - exactly 3 Ts
      - preceded by CCC
      - remove 1 T
    """
    out = []
    i = 0
    n = len(seq)

    while i < n:
        if seq[i].upper() == "T":
            # find full T run
            j = i
            while j < n and seq[j].upper() == "T":
                j += 1

            run_len = j - i
            remove = 0

            # ---------- Type 3 ----------
            if run_len >= 5:
                # index of base preceding the *rightmost* 5 Ts
                idx = i + run_len - 5 - 1
                if idx >= 0 and seq[idx].upper() in "ACT":
                    remove = 2

            # ---------- Type 2 ----------
            elif run_len == 4:
                if i > 0 and seq[i - 1].upper() == "G":
                    remove = 1

            # ---------- Type 1 ----------
            elif run_len == 3:
                if i >= 3 and seq[i - 3:i].upper() == "CCC":
                    remove = 1

            keep = run_len - remove
            if keep > 0:
                out.extend(seq[i:i + keep])

            i = j
        else:
            out.append(seq[i])
            i += 1

    return "".join(out)

# ---------------------------------------------------------
# FASTA I/O
# ---------------------------------------------------------
with open(args.input) as fin, open(args.output, "w") as fout:
    header = None
    seq_buf = []

    for line in fin:
        line = line.rstrip()
        if line.startswith(">"):
            if header:
                unedited = unedit_sequence("".join(seq_buf))
                fout.write(header + "\n")
                if args.wrap:
                    for l in wrap(unedited, args.wrap):
                        fout.write(l + "\n")
                else:
                    fout.write(unedited + "\n")
            header = line
            seq_buf = []
        else:
            seq_buf.append(line)

    # last record
    if header:
        unedited = unedit_sequence("".join(seq_buf))
        fout.write(header + "\n")
        if args.wrap:
            for l in wrap(unedited, args.wrap):
                fout.write(l + "\n")
        else:
            fout.write(unedited + "\n")

print("Unediting complete.")

