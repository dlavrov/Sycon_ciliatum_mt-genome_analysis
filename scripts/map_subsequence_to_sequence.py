#!/usr/bin/env python3

import sys
from Bio import SeqIO
from Bio.Seq import Seq

def find_all(haystack, needle):
    """Yield all start positions of needle in haystack (0-based)."""
    start = 0
    while True:
        pos = haystack.find(needle, start)
        if pos == -1:
            return
        yield pos
        start = pos + 1

def main(orfs_fasta, genome_fasta, out_gff):
    orfs = list(SeqIO.parse(orfs_fasta, "fasta"))
    genome = list(SeqIO.parse(genome_fasta, "fasta"))

    with open(out_gff, "w") as out:
        out.write("##gff-version 3\n")

        for ref in genome:
            ref_seq = str(ref.seq).upper()
            ref_len = len(ref_seq)

            for orf in orfs:
                orf_seq = str(orf.seq).upper()
                orf_len = len(orf_seq)

                # ---------- forward strand ----------
                for pos in find_all(ref_seq, orf_seq):
                    start = pos + 1
                    end = start + orf_len - 1
                    attrs = f"ID={orf.id};Name={orf.id}"
                    out.write(
                        f"{ref.id}\tORF_mapper\tCDS\t"
                        f"{start}\t{end}\t.\t+\t0\t{attrs}\n"
                    )

                # ---------- reverse strand ----------
                rc_seq = str(Seq(orf_seq).reverse_complement())
                for pos in find_all(ref_seq, rc_seq):
                    start = pos + 1
                    end = start + orf_len - 1
                    attrs = f"ID={orf.id}_rc;Name={orf.id}"
                    out.write(
                        f"{ref.id}\tORF_mapper\tCDS\t"
                        f"{start}\t{end}\t.\t-\t0\t{attrs}\n"
                    )

    print(f"GFF3 written to {out_gff}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(
            "Usage: map_orfs_to_gff.py orfs.fasta genome.fasta output.gff"
        )

    main(sys.argv[1], sys.argv[2], sys.argv[3])
