#!/usr/bin/env python3

import re
import argparse
import natsort  # pip install natsort


def parse_seq_id(line: str) -> str:
    """Extract seq_id from '### Results for sequence:' line."""
    seq_id = line.split(":", 1)[1].strip()
    seq_id = seq_id.split()[0]
    seq_id = seq_id.rstrip(";")
    return seq_id


def convert_etandem_to_sorted_gff3(
    input_file: str,
    output_file: str,
    units: bool = False,
    tolerant: bool = False,
    cluster_type: str = "tandem_repeat",
    unit_type: str = "tandem_repeat_unit",
) -> None:
    """
    Convert EMBOSS etandem output to sorted GFF3.

    Default behavior is STRICT:
      require (end-start+1) == (Size * Copies)

    If tolerant=True:
      - in unit mode, emit as many full units as fit
      - in cluster mode, always emit the cluster feature
    """
    gff_records = []

    with open(input_file) as infile:
        seq_id = None
        in_data_block = False

        for line in infile:
            line = line.rstrip()

            if line.startswith("### Results for sequence:"):
                seq_id = parse_seq_id(line)
                continue

            if re.match(r"\s*Start\s+End\s+Strand\s+Score", line):
                in_data_block = True
                continue

            if in_data_block and (line.startswith("#") or line.strip() == ""):
                in_data_block = False
                continue

            if not in_data_block:
                continue

            if not re.match(r"^\s*\d+\s+\d+\s+[+-]", line):
                continue

            parts = re.split(r"\s+", line.strip(), maxsplit=7)
            if len(parts) < 8:
                continue

            try:
                start = int(parts[0])
                end = int(parts[1])
                strand = parts[2]
                score = parts[3]
                size = int(parts[4])
                copies = int(parts[5])
                identity = parts[6]
                consensus = parts[7].replace(";", "%3B")
            except ValueError:
                continue

            if seq_id is None:
                continue

            array_len = end - start + 1
            expected_len = size * copies

            if not tolerant and array_len != expected_len:
                import sys
                print(
                    f"WARNING: skipping {seq_id}:{start}-{end} "
                    f"(span={array_len} != size*copies={expected_len}; "
                    f"size={size}; copies={copies})",
                    file=sys.stderr
                )
                continue

            if not units:
                # Cluster-level feature
                record_line = (
                    f"{seq_id}\tetandem\t{cluster_type}\t{start}\t{end}\t{score}\t{strand}\t.\t"
                    f"ID=tr_{seq_id}_{start}_{end};"
                    f"Size={size};Copies={copies};Identity={identity};Consensus={consensus}\n"
                )
                gff_records.append((seq_id, start, record_line))
            else:
                # Unit-level features
                parent_id = f"tr_{seq_id}_{start}_{end}"

                if tolerant:
                    max_units = min(copies, array_len // size)
                else:
                    max_units = copies

                for i in range(max_units):
                    unit_start = start + i * size
                    unit_end = unit_start + size - 1

                    if unit_end > end:
                        break

                    record_line = (
                        f"{seq_id}\tetandem\t{unit_type}\t{unit_start}\t{unit_end}\t{score}\t{strand}\t.\t"
                        f"ID={parent_id}_copy{i+1};"
                        f"Parent={parent_id};"
                        f"Copy={i+1};"
                        f"Size={size};"
                        f"Copies={copies};"
                        f"Identity={identity};"
                        f"Consensus={consensus}\n"
                    )
                    gff_records.append((seq_id, unit_start, record_line))

    gff_records_sorted = natsort.natsorted(gff_records, key=lambda x: (x[0], x[1]))

    with open(output_file, "w") as out:
        out.write("##gff-version 3\n")
        for _, _, line in gff_records_sorted:
            out.write(line)


def main():
    parser = argparse.ArgumentParser(
        description="Convert EMBOSS etandem output to sorted GFF3 (clusters or split repeat units)."
    )
    parser.add_argument("input_file", help="Path to etandem output file")
    parser.add_argument("output_file", help="Path to output GFF3 file")
    parser.add_argument(
        "--units",
        action="store_true",
        help="Output individual repeat units (one feature per copy) instead of clusters"
    )
    parser.add_argument(
        "--tolerant",
        action="store_true",
        help="Allow span != Size*Copies (emit partial units if needed)"
    )
    parser.add_argument(
        "--cluster-type",
        default="tandem_repeat",
        help="GFF3 type for cluster features (default: tandem_repeat)"
    )
    parser.add_argument(
        "--unit-type",
        default="tandem_repeat_unit",
        help="GFF3 type for unit features (default: tandem_repeat_unit)"
    )

    args = parser.parse_args()

    convert_etandem_to_sorted_gff3(
        args.input_file,
        args.output_file,
        units=args.units,
        tolerant=args.tolerant,
        cluster_type=args.cluster_type,
        unit_type=args.unit_type,
    )


if __name__ == "__main__":
    main()

