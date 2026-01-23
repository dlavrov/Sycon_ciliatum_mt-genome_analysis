#!/usr/bin/env python3
import sys
import re
from collections import defaultdict

def usage():
    print(
        "Usage:\n"
        "  update_gff_strand_from_cluster_blast.py <input.gff> <blast.tsv> <output.gff>\n\n"
        "Assumptions:\n"
        "  - input.gff contains tandem_repeat features with Parent=tr_<seqid>_<start>_<end>\n"
        "  - BLAST TSV is outfmt 6 and includes 'sstrand' as column 7 (0-based index 6)\n"
        "  - BLAST sseqid is like: Contig1:247-667 (i.e., the cluster interval)\n",
        file=sys.stderr
    )
    sys.exit(1)

def parse_cluster_family_from_sseqid(sseqid):
    """
    Convert BLAST subject id like 'Contig1:247-667' into GFF family id:
      'tr_Contig1_247_667'
    Also tolerates extra decorations after spaces or '::'.
    """
    sseqid = sseqid.split("::", 1)[0]
    sseqid = sseqid.split()[0]

    m = re.match(r"^([^:]+):(\d+)-(\d+)$", sseqid)
    if not m:
        return None
    contig, start, end = m.group(1), m.group(2), m.group(3)
    return f"tr_{contig}_{start}_{end}"

def parse_blast_cluster_level(blast_file):
    """
    Returns:
      fam_has_hit: set(family_id)
      fam_all_minus: dict(family_id -> bool)
    Rule:
      - fam_has_hit if any hit exists
      - fam_all_minus is True only if >=1 hit exists and ALL hits are minus
    """
    fam_strands = defaultdict(list)

    with open(blast_file) as bf:
        for line in bf:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 7:
                continue

            sseqid = fields[1]
            sstrand = fields[6].lower()   # 'plus' / 'minus'
            strand = "-" if sstrand == "minus" else "+"

            fam = parse_cluster_family_from_sseqid(sseqid)
            if fam is None:
                continue

            fam_strands[fam].append(strand)

    fam_has_hit = set(fam_strands.keys())
    fam_all_minus = {}

    for fam, strands in fam_strands.items():
        fam_all_minus[fam] = (len(strands) > 0 and all(s == "-" for s in strands))

    return fam_has_hit, fam_all_minus

def update_gff(gff_in, fam_has_hit, fam_all_minus, gff_out):
    """
    - If feature belongs to family with fam_all_minus=True -> set strand to '-'
    - If family has any hit -> add BlastHit=1 attribute (optional but useful)
    """
    with open(gff_in) as infile, open(gff_out, "w") as out:
        for line in infile:
            if line.startswith("#") or not line.strip():
                out.write(line)
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                out.write(line)
                continue

            attrs = parts[8]

            # family id is Parent=... for your individual repeat features
            m = re.search(r"(?:^|;)Parent=([^;]+)", attrs)
            fam = m.group(1) if m else None

            if fam and fam in fam_all_minus and fam_all_minus[fam]:
                parts[6] = "-"

            # tag hit families (helps for plotting color)
            if fam and fam in fam_has_hit:
                if "BlastHit=" not in attrs:
                    parts[8] = attrs + ";BlastHit=1"
                else:
                    parts[8] = attrs

            out.write("\t".join(parts) + "\n")

def main():
    if len(sys.argv) != 4:
        usage()

    gff_in = sys.argv[1]
    blast_tsv = sys.argv[2]
    gff_out = sys.argv[3]

    fam_has_hit, fam_all_minus = parse_blast_cluster_level(blast_tsv)
    update_gff(gff_in, fam_has_hit, fam_all_minus, gff_out)

if __name__ == "__main__":
    main()

