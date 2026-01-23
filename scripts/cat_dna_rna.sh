#!/usr/bin/env bash

for gene in atp6 atp9 cob cox1 cox2 cox3 nad1 nad2 nad3 nad4 nad5; do
    cat "./published/${gene}_published.fa" "./edited_genes/${gene}.lst.fa" "../../data/transcriptome/rnaSpades_trimmomatic/${gene}_transcript.lst.out" "../../data/transcriptome/Trinity_trimmomatic/${gene}_transcript.lst.out" > "./rna_dna_alignments/${gene}_dna_rna.fa"
done
