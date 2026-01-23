#!/usr/bin/env python3

import sys
import re
from collections import defaultdict

def usage():
    print(
        "Usage:\n"
        "  update_gff_strand_from_blast.py <input.gff> <blast.tsv> <output.gff>\n\n"
        "BLAST file must be outfmt 6 and include the 'sstrand' column.\n"
        "The sseqid must contain the GFF feature ID (e.g. from bedtools -name).",
        file=sys.stderr
    )
    sys.exit(1)


def parse_blast(blast_file):
    """
    Returns:
      copy_strand: dict {copy_id -> '+'|'-'}
      family_strand_votes: dict {family_id -> set(['+','-'])}
    """
    copy_strand = {}
    family_votes = defaultdict(set)

    with open(blast_file) as bf:
        for line in bf:
            if not line.strip():
                continue

            fields = line.rstrip().split("\t")
            if len(fields) < 7:
                continue

            sseqid = fields[1]
            sstrand = fields[6]

            strand = "-" if sstrand == "minus" else "+"

            # Extract copy ID before :: if present
            copy_id = sseqid.split("::")[0]

            # Infer family ID by removing _copyN
            family_id = re.sub(r"_copy\d+$", "", copy_id)

            copy_strand[copy_id] = strand
            family_votes[family_id].add(strand)

    return copy_strand, family_votes


def decide_family_strand(family_votes):
    """
    Rule:
      - If any minus hit exists → '-'
      - Else if only plus hits → '+'
    """
    family_strand = {}

    for fam, strands in family_votes.items():
        if "-" in strands:
            family_strand[fam] = "-"
        else:
            family_strand[fam] = "+"

    return family_strand


def update_gff(gff_file, family_strand, output_file):
    with open(gff_file) as infile, open(output_file, "w") as out:
        for line in infile:
            if line.startswith("#"):
                out.write(line)
                continue

            parts = line.rstrip().split("\t")
            if len(parts) != 9:
                out.write(line)
                continue

            attrs = parts[8]

            # Extract Parent= or ID=
            parent_match = re.search(r"Parent=([^;]+)", attrs)
            id_match = re.search(r"ID=([^;]+)", attrs)

            family_id = None

            if parent_match:
                family_id = parent_match.group(1)
            elif id_match:
                family_id = re.sub(r"_copy\d+$", "", id_match.group(1))

            if family_id and family_id in family_strand:
                parts[6] = family_strand[family_id]

            out.write("\t".join(parts) + "\n")


def main():
    if len(sys.argv) != 4:
        usage()

    gff_file = sys.argv[1]
    blast_file = sys.argv[2]
    output_file = sys.argv[3]

    copy_strand, family_votes = parse_blast(blast_file)
    family_strand = decide_family_strand(family_votes)
    update_gff(gff_file, family_strand, output_file)


if __name__ == "__main__":
    main()

