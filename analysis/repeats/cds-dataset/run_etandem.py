#!/usr/bin/env python3
import os
import re
import subprocess
import tempfile
from pathlib import Path
import natsort
import argparse

# ---------- PART 1: Run etandem per sequence ----------
def run_etandem_per_sequence(input_fasta, etandem_output, minrepeat=10, maxrepeat=200):
    """
    Runs etandem on each sequence in a multi-FASTA file and appends results to a combined output.
    """
    input_fasta = Path(input_fasta)
    etandem_output = Path(etandem_output)

    # Read all sequences from the multi-FASTA
    sequences = {}
    seq_id = None
    seq_lines = []
    with open(input_fasta) as infile:
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:
                    sequences[seq_id] = "".join(seq_lines)
                seq_id = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if seq_id:
            sequences[seq_id] = "".join(seq_lines)

    print(f"Found {len(sequences)} sequences in {input_fasta.name}")

    # Clear output file
    etandem_output.write_text("")

    # Temporary file for single sequences
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fa") as temp_fa:
        temp_fasta_path = Path(temp_fa.name)

    # Process each sequence
    for seq_id, seq in sequences.items():
        print(f"Running etandem on sequence: {seq_id}")

        # Write to temp fasta
        with open(temp_fasta_path, "w") as tempf:
            tempf.write(f">{seq_id}\n{seq}\n")

        # Run etandem via subprocess
        try:
            result = subprocess.run(
                [
                    "etandem",
                    "-sequence", str(temp_fasta_path),
                    "-outfile", "stdout",
                    "-minrepeat", str(minrepeat),
                    "-maxrepeat", str(maxrepeat)
                ],
                capture_output=True,
                text=True,
                check=True
            )
        except subprocess.CalledProcessError as e:
            print(f"Error running etandem on {seq_id}: {e}")
            continue

        # Append results
        with open(etandem_output, "a") as out:
            out.write(f"### Results for sequence: {seq_id} ###\n")
            out.write(result.stdout)
            out.write("\n\n")

    os.remove(temp_fasta_path)
    print(f"All etandem results written to {etandem_output}")


# ---------- PART 2: Convert etandem output to GFF3 ----------
def convert_etandem_to_sorted_gff3(input_file, output_file):
    with open(input_file) as infile:
        seq_id = None
        in_data_block = False
        gff_records = []

        for line in infile:
            line = line.rstrip()

            if line.startswith("### Results for sequence:"):
                seq_id = line.split(":")[1].strip()
                continue

            if re.match(r"\s*Start\s+End\s+Strand\s+Score", line):
                in_data_block = True
                continue

            if in_data_block and (line.startswith("#") or line.strip() == ""):
                in_data_block = False
                continue

            if in_data_block:
                if not re.match(r'^\s*\d+\s+\d+\s+[+-]', line):
                    continue

                parts = re.split(r'\s+', line.strip(), maxsplit=7)
                if len(parts) < 8:
                    continue

                try:
                    start = int(parts[0])
                    end = int(parts[1])
                    strand = parts[2]
                    score = parts[3]
                    size = parts[4]
                    count = parts[5]
                    identity = parts[6]
                    consensus = parts[7].replace(";", "%3B")
                except ValueError:
                    continue

                record = (
                    seq_id,
                    start,
                    f"{seq_id}\tetandem\ttandem_repeat\t{start}\t{end}\t{score}\t{strand}\t.\t"
                    f"Note=Size:{size};Copies:{count};Identity:{identity};Consensus={consensus}\n"
                )
                gff_records.append(record)

    # Natural sort by sequence ID, then start coordinate
    gff_records_sorted = natsort.natsorted(gff_records, key=lambda x: (x[0], x[1]))

    with open(output_file, "w") as out:
        out.write("##gff-version 3\n")
        for _, _, line in gff_records_sorted:
            out.write(line)

    print(f"GFF3 file written to {output_file}")


# ---------- MAIN ----------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run etandem on all sequences and convert results to GFF3.")
    parser.add_argument("input_fasta", help="Path to multi-FASTA input file")
    parser.add_argument("etandem_output", help="Path for intermediate etandem results")
    parser.add_argument("final_gff", help="Path for final GFF3 output")
    parser.add_argument("--minrepeat", type=int, default=10, help="Minimum repeat size for etandem")
    parser.add_argument("--maxrepeat", type=int, default=500, help="Maximum repeat size for etandem")
    args = parser.parse_args()

    run_etandem_per_sequence(args.input_fasta, args.etandem_output, args.minrepeat, args.maxrepeat)
    convert_etandem_to_sorted_gff3(args.etandem_output, args.final_gff)

    print("Pipeline complete!")

