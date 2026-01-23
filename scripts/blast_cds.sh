#!/bin/bash
set -ex #exit as soon as any line fails and print each command
query="published_cds_no_editing.fa"
outfile="reads_with_genes.out"
module load blast-plus/2.13.0-py310-irl7uxo
#makeblastdb -in /ptmp/LAS/dlavrov-lab/ERR12319363.fasta -out ERR12319363_blast_db -dbtype nucl -parse_seqids
blastn  -query ${query} -db ERR12319363_blast_db -out ${outfile} -num_threads 16 -evalue 1e-10 -outfmt  "6 qseqid sseqid sskingdoms staxids sscinames scomnames qcovs qcovhsp qcovus evalue bitscore"

