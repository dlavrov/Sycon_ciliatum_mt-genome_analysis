#!/bin/bash
    set -e
    set -u
    set -o pipefail
    
## create databases
mmseqs createdb cc_oc_proteins.aa proteins
mmseqs createdb ../../cap-assembly/all_hits_cap_edited_strand_plus.fa genome_plus
mmseqs createdb ../../cap-assembly/all_hits_cap_edited_strand_minus.fa genome_minus
mmseqs createdb ~/work/projects/Sycon_ciliatum_analysis/data/transcriptome/ERR12318596_spades_rna/transcripts.fasta  transcripts 
## run searches
mmseqs search --translation-table 4 --start-sens 2 -s 7 --sens-steps 3 -a 1 --num-iterations 2 proteins genome_minus proteins_vs_genome_minus tmp
mmseqs search --translation-table 4 --start-sens 2 -s 7 --sens-steps 3 -a 1 --num-iterations 2 proteins genome_plus proteins_vs_genome_plus tmp
mmseqs search --translation-table 4 --start-sens 2 -s 7 --sens-steps 3 -a 1 --num-iterations 2 proteins transcripts proteins_vs_transcripts tmp

## format results
mmseqs convertalis proteins genome_plus proteins_vs_genome_plus proteins_vs_genome_plus.tsv --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits,cigar,qheader,theader,qaln,taln --search-type 2
mmseqs convertalis proteins genome_minus proteins_vs_genome_minus proteins_vs_genome_minus.tsv --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits,cigar,qheader,theader,qaln,taln --search-type 2
mmseqs convertalis proteins transcripts proteins_vs_transcripts proteins_vs_transcripts.tsv --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits,cigar,qheader,theader,qaln,taln --search-type 2
