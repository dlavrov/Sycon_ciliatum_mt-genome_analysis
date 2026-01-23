# README
## Mapping cds
Note: the script I use here maps exact sequence which I took from mf file
```
python map_subsequence_to_sequence.py trinity_best_matched_cds.fa trinity_best_matches_in_genome.fa cds.gff
```
## Mapping etandem repeats
This script extracts individual etandem repeats from etandem outfile
The old version would provide coordinates for the whole repeat unit
```
~/work/projects/Sycon_ciliatum_analysis/Sycon_ciliatum_mt-genome_final/scripts/etandem2gff3.py
etandem2gff3.py trinity_best_matches_in_genome_etandem.out trinity_best_matches_in_genome_etandem.gff
etandem2gff3_cluster.py trinity_best_matches_in_genome_etandem.out trinity_best_matches_in_genome_etandem_runs.gff
```

## Extracting repeats from chromosomes
```
module load bedtools2/2.31.1-py311-6kemgt3
bedtools getfasta -fi trinity_best_matches_in_genome.fa -bed trinity_best_matches_in_genome_etandem.gff -fo tandem_repeats.fa
bedtools getfasta -fi trinity_best_matches_in_genome.fa -bed trinity_best_matches_in_genome_etandem_runs.gff -fo tandem_repeats_runs.fa
```

## Trying to put strand orientation on the repeats
```
module load blast-plus/2.13.0-py310-irl7uxo
# makeblastdb -in tandem_repeats.fa -dbtype nucl -out repeat_db
# makeblastdb -in tandem_repeats_runs.fa -dbtype nucl -out repeat_runs.db
blastn -query known_repeats.fa -db repeat_db -out repeats_vs_repeats.tsv -outfmt "6 qseqid sseqid pident length qlen slen sstrand qstart qend sstart send bitscore" -evalue 1e-5
./update_gff_strand_from_blast.py trinity_best_matches_in_genome_etandem.gff repeats_vs_repeats.tsv test.gff
blastn -query known_repeats.fa -db repeat_runs.db -out repeats_vs_repeat_runs.tsv -outfmt "6 qseqid sseqid pident length qlen slen sstrand qstart qend sstart send bitscore" -evalue 1e-5
./update_gff_strand_from_blast.py trinity_best_matches_in_genome_etandem.gff repeats_vs_repeat_runs.tsv repeats_in_runs.gff
```
## Coloring repeats
```
./polyAT_to_gff.py trinity_best_matches_in_genome.fa AT.gff --min-len 10 --end-window 100
cat AT.gff >> test.gff
./update_gff_strand_from_blast_new.py etandem_repeats.gff repeats_vs_repeat_runs.tsv repeats_in_runs.gff

module load py-matplotlib/3.10.1-py311-luy2jpy
module load py-pandas/2.2.3-py311-zlbcman
module load py-biopython/1.81-py311-xigmold

./gff2pdf_final.py --blast repeats_vs_repeats.tsv --align left test.gff trinity_best_matches_in_genome.fa chromosomes_best_Trinity_matches.pdf
./gff2pdf_final_updated.py --blast repeats_vs_repeat_runs.tsv test.gff trinity_best_matches_in_genome.fa test.pdf
```

Here is a history of updates to my mapping script
- Dec 30 14:20 gff2pdf_final_updated.py
- Dec 30 13:30 update_gff_strand_from_blast_new.py
- Dec 24 14:39 irf2gff.py
- Dec 24 10:52 gff2pdf_final.py
- Dec 24 10:50 gff2pdf_save.py
- Dec 24 10:09 gff2pdf_cvg.py
- Dec 23 22:21 gff2pdf_test.py
- Dec 23 22:18 gff2pdf_pub.py
- Dec 23 22:01 polyAT_to_gff.py
- Dec 23 21:30 gff2pdf_new.py
- Dec 23 21:11 gff2pdf.py
- Dec 23 20:40 update_gff_strand_from_blast.py
- Dec 23 17:55 merge_blast_gff.py
- Dec 23 17:34 blast_to_gff.py
- Dec 22 19:22 gff2pdf_old.py
- Dec 22 18:38 map_subsequence_to_sequence.py
- Dec 22 18:32 de-edit_transcripts.py
- Dec 22 14:58 gff2pdf_works.py
