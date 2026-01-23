#!/usr/bin/env python3

import argparse

def blast_to_gff(blast_file, output_gff, feature_type="blast_hit"):
    with open(blast_file) as fin, open(output_gff, "w") as fout:
        fout.write("##gff-version 3\n")

        hit_id = 0
        for line in fin:
            if not line.strip():
                continue

            (
                qseqid, sseqid, pident, length, mismatch, gapopen,
                qstart, qend, sstart, send, evalue, bitscore
            ) = line.strip().split("\t")

            sstart = int(sstart)
            send = int(send)
            qstart = int(qstart)
            qend = int(qend)

            if sstart <= send:
                start, end = sstart, send
                strand = "+"
            else:
                start, end = send, sstart
                strand = "-"

            hit_id += 1
            attrs = (
                f"ID=blast_{hit_id};"
                f"Target={qseqid} {qstart} {qend};"
                f"Identity={pident};"
                f"Length={length};"
                f"Evalue={evalue}"
            )

            fout.write(
                f"{sseqid}\tblast\t{feature_type}\t"
                f"{start}\t{end}\t{bitscore}\t"
                f"{strand}\t.\t{attrs}\n"
            )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert BLAST tabular output (outfmt 6) to GFF3"
    )
    parser.add_argument("blast", help="BLAST output (outfmt 6)")
    parser.add_argument("gff", help="Output GFF3 file")
    parser.add_argument(
        "--type",
        default="blast_hit",
        help="GFF feature type (default: blast_hit)"
    )
    args = parser.parse_args()

    blast_to_gff(args.blast, args.gff, args.type)

