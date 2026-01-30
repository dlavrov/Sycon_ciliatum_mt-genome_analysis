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
parser.add_argument(
    "-o", "--output",
    help="Output FASTA/masterfile (default: stdout)"
)
parser.add_argument(
    "--mode",
    choices=["inline", "master"],
    default="inline",
    help="Output mode: 'inline' = edited FASTA only; 'master' = edited sequence with interleaved comments."
)
parser.add_argument(
    "--gff",
    help="If provided, write GFF3 annotations to this file (genomic coordinates)."
)
parser.add_argument(
    "--wrap", type=int, default=0,
    help="FASTA line width (0 = no wrapping). In master mode, wrapping applies to sequence chunks between comment blocks."
)
args = parser.parse_args()
want_gff = bool(args.gff)

# ---------------------------------------------------------
# Annotation helpers
# ---------------------------------------------------------
def master_comment_block(type_no: int, idx: int, motif_no_insertions: str):
    tag = f"G-Mot-editing_site_type{type_no}_{idx}"
    return [
        f";     {tag} ==> start",
        motif_no_insertions,  # motif WITHOUT inserted t's
        f";     {tag} ==> end",
    ]

def gff_entry(seqid, start_g, end_g, type_no, idx, inserted, context, edited_start, edited_end):
    # start_g/end_g are genomic coords (original sequence, 1-based, inclusive)
    # edited_start/edited_end are coords on the edited sequence for the GENOMIC motif only
    # (inserted t's excluded from the span), but shifted by earlier insertions
    attrs = (
        f"ID=editing_site_type{type_no}_{idx};"
        f"Type=type{type_no};"
        f"InsertedTs={inserted};"
        f"EditedStart={edited_start};"
        f"EditedEnd={edited_end};"
        f"Context={context}"
    )
    return (
        f"{seqid}\tmt_editing\tediting_site_type{type_no}\t"
        f"{start_g}\t{end_g}\t.\t+\t.\t{attrs}"
    )

def write_seq_chunk(lines_out, chunk, wrap_width):
    if not chunk:
        return
    if wrap_width and wrap_width > 0:
        lines_out.extend(wrap(chunk, wrap_width))
    else:
        lines_out.append(chunk)

# ---------------------------------------------------------
# Core processing
# ---------------------------------------------------------
def process_sequence(seq, seqid, wrap_width, mode, gff_out, want_gff):
    out_lines = []
    counters = {1: 0, 2: 0, 3: 0}
    inserted_total = 0

    # Total inserted bases BEFORE current genomic position
    ins_before = 0

    # Output buffers
    inline_parts = []     # inline mode: build edited sequence
    chunk_buf = []        # master mode: sequence between comment blocks

    def flush_master_chunk():
        nonlocal chunk_buf
        if chunk_buf:
            write_seq_chunk(out_lines, "".join(chunk_buf), wrap_width)
            chunk_buf = []

    i = 0
    while i < len(seq):

        if seq[i] == "T":
            j = i
            while j < len(seq) and seq[j] == "T":
                j += 1

            run_len = j - i
            type_no = None
            insert = 0
            prefix_len = 0

            # ≥4 Ts → type 3 (insert 2)
            if run_len >= 4:
                type_no, insert, prefix_len = 3, 2, 0

            # exactly 3 Ts
            elif run_len == 3 and i > 0:
                prev = seq[i - 1]
                if prev in "ACT":
                    type_no, insert, prefix_len = 3, 2, 1  # {A/C/T}TTT
                elif prev == "G":
                    type_no, insert, prefix_len = 2, 1, 1  # GTTT

            # exactly 2 Ts with CCC upstream
            elif run_len == 2 and i >= 3 and seq[i - 3:i] == "CCC":
                type_no, insert, prefix_len = 1, 1, 3      # CCCTT

            if type_no:
                counters[type_no] += 1
                inserted_total += insert

                motif_prefix = seq[i - prefix_len:i] if prefix_len else ""
                motif_no_insertions = motif_prefix + ("T" * run_len)

                # Genomic coords include full motif (prefix + T-run), NOT inserted t's
                start_g = (i - prefix_len) + 1
                end_g = i + run_len
                motif_genomic_len = end_g - start_g + 1

                # Edited coords span ONLY the genomic motif (exclude inserted t's),
                # shifted by prior insertions
                edited_start = start_g + ins_before
                edited_end = edited_start + motif_genomic_len - 1

                if want_gff:
                    context = motif_no_insertions  # motif only (no inserted t's)
                    gff_out.append(
                        gff_entry(
                            seqid, start_g, end_g,
                            type_no, counters[type_no],
                            insert, context,
                            edited_start, edited_end
                        )
                    )

                if mode == "master":
                    # Move prefix bases into the motif block: remove them from upstream chunk
                    if prefix_len and len(chunk_buf) >= prefix_len:
                        del chunk_buf[-prefix_len:]
                    elif prefix_len:
                        chunk_buf = []

                    flush_master_chunk()

                    # Comment block contains motif (no insertions) as the motif line
                    out_lines.extend(
                        master_comment_block(type_no, counters[type_no], motif_no_insertions)
                    )

                    # After the motif block, output ONLY the inserted t's
                    if insert > 0:
                        out_lines.append("t" * insert)

                else:
                    # Inline mode: prefix already emitted earlier in stream,
                    # so add only T-run + inserted t's
                    inline_parts.append(("T" * run_len) + ("t" * insert))

                # Downstream edited coords shift by 'insert'
                ins_before += insert

                i = j
                continue

            # Not an editing site: emit one T
            if mode == "master":
                chunk_buf.append(seq[i])
            else:
                inline_parts.append(seq[i])
            i += 1
            continue

        # non-T base
        if mode == "master":
            chunk_buf.append(seq[i])
        else:
            inline_parts.append(seq[i])
        i += 1

    # flush remaining
    if mode == "master":
        flush_master_chunk()
        edited_len = len(seq) + inserted_total
    else:
        full_edited = "".join(inline_parts)
        edited_len = len(full_edited)
        out_lines = wrap(full_edited, wrap_width) if wrap_width else [full_edited]

    summary = {
        "seqid": seqid,
        "type1": counters[1],
        "type2": counters[2],
        "type3": counters[3],
        "inserted_total": inserted_total,
        "orig_len": len(seq),
        "edited_len": edited_len,
    }
    return out_lines, summary

# ---------------------------------------------------------
# Summary helpers (stderr)
# ---------------------------------------------------------
def write_summary_header():
    print(
        "seqid\torig_len\tedited_len\ttype1\ttype2\ttype3\tinserted_t_total",
        file=sys.stderr
    )

def write_summary_row(s):
    print(
        f"{s['seqid']}\t{s['orig_len']}\t{s['edited_len']}\t"
        f"{s['type1']}\t{s['type2']}\t{s['type3']}\t{s['inserted_total']}",
        file=sys.stderr
    )

# ---------------------------------------------------------
# FASTA I/O
# ---------------------------------------------------------
gff_lines = ["##gff-version 3"]
summaries = []
totals = {"type1": 0, "type2": 0, "type3": 0, "inserted_total": 0}

out_handle = open(args.output, "w") if args.output else sys.stdout

try:
    with open(args.input) as fin:
        header = None
        seq_parts = []

        for line in fin:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header:
                    seqid = header[1:].split()[0]
                    out_lines, summary = process_sequence(
                        "".join(seq_parts), seqid, args.wrap, args.mode, gff_lines, want_gff
                    )

                    out_handle.write(header + "\n")
                    for l in out_lines:
                        out_handle.write(l + "\n")

                    summaries.append(summary)
                    totals["type1"] += summary["type1"]
                    totals["type2"] += summary["type2"]
                    totals["type3"] += summary["type3"]
                    totals["inserted_total"] += summary["inserted_total"]

                header = line
                seq_parts = []
            else:
                if line and not line.startswith(";"):
                    seq_parts.append(line.strip())

        # last record
        if header:
            seqid = header[1:].split()[0]
            out_lines, summary = process_sequence(
                "".join(seq_parts), seqid, args.wrap, args.mode, gff_lines, want_gff
            )

            out_handle.write(header + "\n")
            for l in out_lines:
                out_handle.write(l + "\n")

            summaries.append(summary)
            totals["type1"] += summary["type1"]
            totals["type2"] += summary["type2"]
            totals["type3"] += summary["type3"]
            totals["inserted_total"] += summary["inserted_total"]

finally:
    if args.output:
        out_handle.close()

# ---------------------------------------------------------
# Write GFF if requested
# ---------------------------------------------------------
if want_gff:
    with open(args.gff, "w") as gff:
        gff.write("\n".join(gff_lines) + "\n")

# ---------------------------------------------------------
# Summary to stderr
# ---------------------------------------------------------
write_summary_header()
for s in summaries:
    write_summary_row(s)

print(
    f"TOTAL\t.\t.\t{totals['type1']}\t{totals['type2']}\t{totals['type3']}\t{totals['inserted_total']}",
    file=sys.stderr
)

