# README
## HMMER search
```
mafft --auto sequences.fasta > aligned.fasta
hmmbuild nad4L.hmm uniprot-nad4L_aln.fa
# Example with TransDecoder
TD2.LongOrfs -G4 -m50 -M50 -t ERR12318596_transcripts.fasta -O ERR12318596_TD2
hmmsearch --tblout results.tbl nad6.hmm ../../data/transcriptome/ERR12318596_TD2/longest_orfs.pep
hmmsearch --tblout results_nad4L.tbl nad4L.hmm ../../data/transcriptome/ERR12318596_TD2/longest_orfs.pep
```

None of the search found any similar transcript.
The best "atp8" seq belongs to the nuclear chromosome 11.

Repeated the search for all mitochondrial chromosomes and found nothing promising.
