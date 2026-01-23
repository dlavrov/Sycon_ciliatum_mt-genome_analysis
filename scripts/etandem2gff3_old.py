#!/usr/bin/env python3

import re
import argparse
import natsort  # pip install natsort

def convert_etandem_to_sorted_gff3(input_file, output_file):
    with open(input_file) as infile:
        seq_id = None
        in_data_block = False
        gff_records = []

        for line in infile:
            line = line.rstrip()

            # Identify sequence/chromosome block
            #if line.startswith("### Results for sequence:"):
            #    seq_id = line.split(":")[1].strip()
            #    continue
            if line.startswith("### Results for sequence:"):
                seq_id = line.split(":", 1)[1].strip()
                seq_id = seq_id.split()[0]      # remove anything after space
                seq_id = seq_id.rstrip(";")     # clean trailing ;;
                continue

            # Detect start of repeat data table
            if re.match(r"\s*Start\s+End\s+Strand\s+Score", line):
                in_data_block = True
                continue

            # Detect end of table or comment line
            if in_data_block and (line.startswith("#") or line.strip() == ""):
                in_data_block = False
                continue

            if in_data_block:
                # Only parse lines starting with numeric start/end and strand (+/-)
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
                    consensus = parts[7].replace(";", "%3B")  # Escape semicolons for GFF
                except ValueError:
                    continue

                record = (
                    seq_id,
                    start,
                    f"{seq_id}\tetandem\ttandem_repeat\t{start}\t{end}\t{score}\t{strand}\t.\t"
                    f"ID=tr_{seq_id}_{start}_{end};"
                    f"Size={size};Copies={count};Identity={identity};Consensus={consensus}\n"
                )
                # record = (
                #    seq_id,
                #    start,
                #    f"{seq_id}\tetandem\ttandem_repeat\t{start}\t{end}\t{score}\t{strand}\t.\t"
                #    f"Note=Size:{size};Copies:{count};Identity:{identity};Consensus={consensus}\n"
                #)
                gff_records.append(record)

    # Sort by chromosome (naturally) and by start position
    gff_records_sorted = natsort.natsorted(gff_records, key=lambda x: (x[0], x[1]))

    # Write GFF3 output
    with open(output_file, 'w') as out:
        out.write("##gff-version 3\n")
        for _, _, line in gff_records_sorted:
            out.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert EMBOSS etandem output to sorted GFF3.")
    parser.add_argument("input_file", help="Path to etandem output file")
    parser.add_argument("output_file", help="Path to output GFF3 file")
    args = parser.parse_args()

    convert_etandem_to_sorted_gff3(args.input_file, args.output_file)

