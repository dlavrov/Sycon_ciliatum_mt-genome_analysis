# README
Analysis of variants found among genes

## Files
- edited_genes: directory with genomic info
	- all80.hits
	- all80_hits_edited.fa: extracted sequences conceptually edited
	- all80_hits.fa: extracted sequences
- fasta_best_hits_published_cds.lst: best hits from fasta search using published mRNA
- genetic_distances: directory to calculate distances among genetic sequences
- pargene_export.out: ml trees based on genomic/transcriptomic sequences
- published: published transcriptomic sequences together and by gene
- README.md: this file
- rna_dna_alignments: alignment created by mafft (original fa sequences deleted)

## Scripts
Script `cat_dna_rna.sh` was used to combine genomic and transcriptomic data from several sources.

## Comments
File gene-contigs-singlets.fa was deleted as it was identical to ../../cap-assembly/mt-cds_cap_all.fa
