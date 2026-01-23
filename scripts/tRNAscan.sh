#!/bin/bash
    set -e
    set -u
    set -o pipefail

infile="$1"
outfile="$(basename "${infile%.*}_tRNAscan")"
tRNAscan-SE -O ${infile} -o ${outfile}.out -f  ${outfile}.str -m  ${outfile}.stat -a  ${outfile}.fa -b  ${outfile}.bed
