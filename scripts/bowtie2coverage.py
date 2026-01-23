#!/usr/bin/env python3

import sys
from collections import defaultdict

def usage():
    sys.exit(f"Usage: {sys.argv[0]} <bowtie1_output.txt>")

if len(sys.argv) != 2:
    usage()

infile = sys.argv[1]

# Nested dictionary:
# data[contig][position][nucleotide] = count
data = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

# Bowtie1 ASCII qualities A–J correspond to low Phred (as in your Perl script)
VALID_QUALS = set("ABCDEFGHIJ")

with open(infile) as fh:
    for line in fh:
        if not line.strip():
            continue

        fields = line.rstrip().split()
        if len(fields) < 5:
            continue

        strand, contig, pos, read, qual = fields[:5]
        pos = int(pos)

        reads = list(read)
        quals = list(qual)

        if strand == "-":
            reads.reverse()
            quals.reverse()

        read_len = len(reads)

        for i, (base, q) in enumerate(zip(reads, quals)):
            if q not in VALID_QUALS:
                continue

            # replicate Perl spos logic exactly
            spos = read_len - 1 - i if strand == "-" else i
            genomic_pos = pos + spos

            data[contig][genomic_pos][base] += 1

# Output
for contig in sorted(data):
    for pos in sorted(data[contig]):
        counts = data[contig][pos]
        line = [
            contig,
            str(pos),
            str(counts.get("A", 0)),
            str(counts.get("G", 0)),
            str(counts.get("C", 0)),
            str(counts.get("T", 0)),
        ]
        print("\t".join(line))

