import argparse
import pandas as pd
import matplotlib.pyplot as plt

# ---------------------------------------------------------
# CLI
# ---------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Plot per-contig coverage from a tab-delimited file"
)
parser.add_argument(
    "coverage",
    help="Coverage file (tab-delimited: contig, position, depth)"
)
parser.add_argument(
    "-o", "--output",
    default="all_contigs_coverage.pdf",
    help="Output PDF file (default: all_contigs_coverage.pdf)"
)
args = parser.parse_args()

# ---------------------------------------------------------
# Load data
# ---------------------------------------------------------
coverage = pd.read_csv(
    args.coverage,
    sep="\t",
    names=["contig", "pos", "depth"]
)

# Ensure numeric and sorted
coverage["pos"] = pd.to_numeric(coverage["pos"])
coverage["depth"] = pd.to_numeric(coverage["depth"])
coverage = coverage.sort_values(["contig", "pos"])

contigs = coverage["contig"].unique()
num_contigs = len(contigs)

# ---------------------------------------------------------
# Plot
# ---------------------------------------------------------
fig, axs = plt.subplots(
    num_contigs,
    1,
    figsize=(12, 4 * num_contigs),
    sharex=False
)

# Ensure axs is iterable
if num_contigs == 1:
    axs = [axs]

for ax, contig in zip(axs, contigs):
    subset = coverage[coverage["contig"] == contig]
    ax.plot(subset["pos"], subset["depth"])
    ax.set_title(f"Coverage for contig: {contig}")
    ax.set_ylabel("Depth")
    ax.grid(True)

axs[-1].set_xlabel("Position")
plt.tight_layout()
plt.savefig(args.output)

print(f"Coverage plot saved as '{args.output}'")

# ---------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------
mean_coverage = coverage.groupby("contig")["depth"].mean()

print("\nMean coverage per contig:")
print(mean_coverage.to_string())

