#!/usr/bin/env python3
import argparse
import re
import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO
from matplotlib.patches import Polygon
from matplotlib.patches import Patch


def add_label_box(ax, x, y, text, fontsize=10):
    ax.text(
        x, y, text,
        ha="center", va="center",
        fontsize=fontsize,
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="black", lw=0.7)
    )


def parse_gff_attrs(attr_str):
    """Parse GFF3 attributes into dict."""
    d = {}
    if not isinstance(attr_str, str):
        return d
    for part in attr_str.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
    return d


def load_blast_matched_repeat_run_keys(blast_tsv):
    """
    Read BLAST outfmt 6 TSV and return a set of subject repeat-run keys.
    Expected subject examples:
      Contig1:247-667
      ERR12319363.1.4671060:5604-6129
    """
    matched = set()
    with open(blast_tsv) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            sseqid = cols[1].split("::", 1)[0]
            if ":" in sseqid and "-" in sseqid:
                matched.add(sseqid)
    return matched


def plot_features(gff_file, fasta_file, outpath, blast_tsv=None,
                  align="left", cds_height=0.9, rep_height=1.0,
                  label_font=13, title_font=14):

    # Read GFF3
    cols = ["seqid","source","type","start","end","score","strand","phase","attrs"]
    df = pd.read_csv(gff_file, sep="\t", comment="#", names=cols)

    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["start","end"])

    # Parse attributes
    df["_attrs"] = df["attrs"].apply(parse_gff_attrs)
    df["ID"] = df["_attrs"].apply(lambda d: d.get("ID", ""))
    df["Parent"] = df["_attrs"].apply(lambda d: d.get("Parent", ""))

    # Load contig lengths
    contig_lengths = {rec.id: len(rec.seq) for rec in SeqIO.parse(fasta_file, "fasta")}
    contigs = list(contig_lengths.keys())
    max_len = max(contig_lengths.values()) if contig_lengths else 0

    # ------------------------------------------------------------
    # Determine which repeat families have BLAST hits (run-based)
    # ------------------------------------------------------------
    hit_families = set()
    if blast_tsv:
        matched_run_keys = load_blast_matched_repeat_run_keys(blast_tsv)
        rep_rows = df[df["type"] == "tandem_repeat"].copy()
        rep_rows["family"] = rep_rows["Parent"].where(
            rep_rows["Parent"] != "", rep_rows["ID"]
        )

        grouped = rep_rows.groupby(["seqid", "family"], dropna=False)
        for (seqid, fam), g in grouped:
            fam_start = int(g["start"].min())
            fam_end = int(g["end"].max())
            run_start0 = fam_start - 1
            run_key = f"{seqid}:{run_start0}-{fam_end}"
            if run_key in matched_run_keys:
                hit_families.add(fam)

    # Colors
    cds_color = "tab:blue"
    repeat_default = "tab:red"
    repeat_hitfam = "tab:green"

    fig, axs = plt.subplots(
        len(contigs), 1,
        figsize=(12, 1.9 * len(contigs)),
        sharex=True
    )
    if len(contigs) == 1:
        axs = [axs]

    for ax, contig in zip(axs, contigs):
        sub = df[df["seqid"] == contig]
        contig_len = contig_lengths[contig]

        offset = (max_len - contig_len) / 2 if align == "center" else 0

        # Backbone
        ax.hlines(0, offset, offset + contig_len, color="black", linewidth=5, zorder=1)

        # Clean axes
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.tick_params(left=False, bottom=False, labelleft=False)
        ax.set_yticks([])

        # CDS
        cds = sub[sub["type"] == "CDS"]
        for _, row in cds.iterrows():
            s = offset + row.start
            e = offset + row.end
            width = e - s

            ax.broken_barh(
                [(s, width)],
                (-cds_height/2, cds_height),
                facecolors=cds_color,
                zorder=2
            )

            if row.ID and width > 120:
                ax.text(
                    s + width / 2, 0, row.ID,
                    ha="center", va="center",
                    fontsize=10, fontstyle="italic",
                    color="white", zorder=4, clip_on=True
                )

        # ------------------------------------------------------------
        # Tandem repeats:
        #   - BLAST-hit families → oriented triangles
        #   - No-hit families   → rhombi (diamonds)
        # ------------------------------------------------------------
        reps = sub[sub["type"] == "tandem_repeat"]
        for _, row in reps.iterrows():
            s = offset + row.start
            e = offset + row.end
            mid = (s + e) / 2
            strand = str(row.strand).strip()
            parent = row.Parent or row.ID

            if parent in hit_families:
                color = repeat_hitfam
                if strand == "-":
                    verts = [(s, 0), (e, rep_height/2), (e, -rep_height/2)]
                else:
                    verts = [(e, 0), (s, rep_height/2), (s, -rep_height/2)]
            else:
                color = repeat_default
                verts = [
                    (mid, rep_height/2),
                    (e, 0),
                    (mid, -rep_height/2),
                    (s, 0)
                ]

            ax.add_patch(
                Polygon(verts, closed=True,
                        facecolor=color, edgecolor="black",
                        linewidth=0.6, zorder=3)
            )

        # polyA / polyT labels
        for _, row in sub[sub["type"] == "polyA_run"].iterrows():
            mid = offset + (row.start + row.end) / 2
            add_label_box(ax, mid, 0.85, "A", fontsize=label_font)

        for _, row in sub[sub["type"] == "polyT_run"].iterrows():
            mid = offset + (row.start + row.end) / 2
            add_label_box(ax, mid, 0.85, "T", fontsize=label_font)

        # Contig label
        if align == "left":
            ax.text(
                offset, 1.05, contig,
                transform=ax.get_xaxis_transform(),
                ha="left", va="bottom",
                fontsize=title_font, fontweight="bold"
            )
        else:
            ax.set_title(contig, fontsize=title_font, fontweight="bold", pad=6)

        ax.set_ylim(-1.15, 1.15)
        ax.set_xlim(0, max_len if max_len else contig_len)

    # Legend
    handles = [
        Patch(facecolor=cds_color, edgecolor="black", label="CDS"),
        Patch(facecolor=repeat_default, edgecolor="black", label="tandem_repeat (no BLAST hit; diamond)"),
    ]
    if blast_tsv:
        handles.append(
            Patch(facecolor=repeat_hitfam, edgecolor="black",
                  label="tandem_repeat (BLAST-hit family; triangle)")
        )

    fig.legend(handles=handles, loc="upper right", frameon=False)

    axs[-1].tick_params(bottom=True, labelbottom=True)
    axs[-1].set_xlabel("Position (bp)", fontsize=12)

    plt.tight_layout()

    if outpath.lower().endswith(".svg"):
        plt.savefig(outpath, format="svg")
    else:
        plt.savefig(outpath)

    print(f"Schematic saved as {outpath}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gff", help="Combined GFF3 (CDS + tandem_repeat + polyA_run/polyT_run)")
    parser.add_argument("fasta", help="Genome FASTA")
    parser.add_argument("out", help="Output figure (.svg for SVG)")
    parser.add_argument("--blast", default=None,
                        help="BLAST outfmt 6 TSV with run subjects (e.g. Contig1:247-667)")
    parser.add_argument("--align", choices=["left", "center"], default="left")

    args = parser.parse_args()

    plot_features(
        args.gff, args.fasta, args.out,
        blast_tsv=args.blast,
        align=args.align
    )

