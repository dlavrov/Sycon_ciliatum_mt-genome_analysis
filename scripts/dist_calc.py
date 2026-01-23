#!/usr/bin/env python3
import sys
import itertools
import statistics
import argparse
from Bio import AlignIO

def calculate_p_distance(seq1, seq2, ignore_gaps=False):
    differences = 0
    compared_sites = 0
    for a, b in zip(seq1, seq2):
        if ignore_gaps and ('-' in (a, b)):
            continue
        if not ignore_gaps and a == '-' and b == '-':
            continue  # skip gap-gap positions only
        compared_sites += 1
        if a != b:
            differences += 1
    if compared_sites == 0:
        return 0.0
    return differences / compared_sites

def compute_distance_matrix(alignment, ignore_gaps=False):
    n = len(alignment)
    matrix = [[0.0]*n for _ in range(n)]
    distances = []
    for i, j in itertools.combinations(range(n), 2):
        p = calculate_p_distance(alignment[i].seq.upper(),
                                 alignment[j].seq.upper(),
                                 ignore_gaps=ignore_gaps)
        matrix[i][j] = matrix[j][i] = p
        distances.append(p)
    return matrix, distances

def print_phylip_matrix(names, matrix, out_stream=sys.stdout):
    print(len(names), file=out_stream)
    for name, row in zip(names, matrix):
        print(f"{name[:10]:<10}", end=' ', file=out_stream)
        for dist in row:
            print(f"{dist:.5f}", end=' ', file=out_stream)
        print(file=out_stream)

def print_stats(distances, err_stream=sys.stderr):
    print("\nBasic Statistics of Pairwise Distances (excluding self):", file=err_stream)
    print(f"Min   : {min(distances):.5f}", file=err_stream)
    print(f"Max   : {max(distances):.5f}", file=err_stream)
    print(f"Mean  : {statistics.mean(distances):.5f}", file=err_stream)
    print(f"Median: {statistics.median(distances):.5f}", file=err_stream)
    if len(distances) > 1:
        print(f"Stdev : {statistics.stdev(distances):.5f}", file=err_stream)
    else:
        print("Stdev : N/A", file=err_stream)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate p-distances from alignment, gap-aware.")
    parser.add_argument("alignment", help="Input alignment file (FASTA)")
    parser.add_argument("--ignore-gaps", action="store_true",
                        help="Ignore any positions containing a gap in either sequence")

    args = parser.parse_args()

    alignment = AlignIO.read(args.alignment, "fasta")
    names = [record.id for record in alignment]

    matrix, distances = compute_distance_matrix(alignment, ignore_gaps=args.ignore_gaps)
    print_phylip_matrix(names, matrix, out_stream=sys.stdout)
    print_stats(distances, err_stream=sys.stderr)

