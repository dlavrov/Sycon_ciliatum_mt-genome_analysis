#!/usr/bin/env python3
import sys
import argparse

def irf_ngs_to_gff3(irf_out, gff_out, source="IRF", write_arms=True):
    print("##gff-version 3", file=gff_out)

    seqid = None
    rep_i = 0

    for raw in irf_out:
        line = raw.strip()
        if not line:
            continue

        # Sequence header lines in -ngs output start with '@'
        # Example: @1 CHROMOSOME ...
        if line.startswith("@"):
            # keep just the first token after '@' if possible; otherwise whole header
            # Many people want the FASTA contig name here; adjust if your @header contains it.
            # In the example page, the header begins with "@1 CHROMOSOME ...".
            # If your IRF preserves FASTA IDs, you might prefer: seqid = line[1:].split()[0]
            seqid = line[1:].split()[0]
            continue

        # Data lines begin with numbers
        parts = line.split()
        if len(parts) < 10:
            continue

        try:
            l_start = int(parts[0])
            l_end   = int(parts[1])
            l_len   = int(parts[2])
            r_start = int(parts[3])
            r_end   = int(parts[4])
            r_len   = int(parts[5])
            loop    = int(parts[6])
            pmatch  = float(parts[7])
            pindel  = float(parts[8])
            score   = float(parts[9])
        except ValueError:
            continue

        if seqid is None:
            # If IRF output didn't include @headers for some reason
            seqid = "unknown"

        rep_i += 1
        rep_id = f"irf_{seqid}_{l_start}_{r_end}_{rep_i}"

        # Full inverted repeat span
        start = min(l_start, l_end, r_start, r_end)
        end   = max(l_start, l_end, r_start, r_end)

        attrs = [
            f"ID={rep_id}",
            f"LeftStart={l_start}",
            f"LeftEnd={l_end}",
            f"LeftLen={l_len}",
            f"RightStart={r_start}",
            f"RightEnd={r_end}",
            f"RightLen={r_len}",
            f"Loop={loop}",
            f"PctMatch={pmatch}",
            f"PctIndel={pindel}",
        ]
        print(
            f"{seqid}\t{source}\tinverted_repeat\t{start}\t{end}\t{score}\t.\t.\t" + ";".join(attrs),
            file=gff_out
        )

        if write_arms:
            # left arm
            print(
                f"{seqid}\t{source}\tIR_arm\t{min(l_start,l_end)}\t{max(l_start,l_end)}\t{score}\t+\t.\t"
                f"ID={rep_id}_L;Parent={rep_id}",
                file=gff_out
            )
            # right arm (reverse-oriented relative to left in an inverted repeat)
            print(
                f"{seqid}\t{source}\tIR_arm\t{min(r_start,r_end)}\t{max(r_start,r_end)}\t{score}\t-\t.\t"
                f"ID={rep_id}_R;Parent={rep_id}",
                file=gff_out
            )

def main():
    ap = argparse.ArgumentParser(
        description="Convert IRF -ngs output (stdout capture) to GFF3."
    )
    ap.add_argument("irf_out", help="IRF output captured with -ngs (e.g., irf.out)")
    ap.add_argument("gff3", help="Output GFF3 file")
    ap.add_argument("--no-arms", action="store_true", help="Do not emit IR_arm features")
    ap.add_argument("--source", default="IRF", help="GFF source field (default: IRF)")
    args = ap.parse_args()

    with open(args.irf_out) as fin, open(args.gff3, "w") as fout:
        irf_ngs_to_gff3(fin, fout, source=args.source, write_arms=(not args.no_arms))

if __name__ == "__main__":
    main()

