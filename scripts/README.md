# README.md
## List of scripts
- annotate\_editing.py: python script to make or annotates predicted editing sites in genomic sequences
- blast.job: 
- blast\_cds.sh: bash script to run command-line blast command
- blast\_to\_gff.py: making gff file for blast hits
- bowtie1.sh: bash script to run bowtie1 command with strict options; calculate
  base-level and mean coverage.
- bowtie1\_only.sh: bash script to only run bowtie1 command
- bowtie2.sh: bash script to run bowtie2 command, calculate per base and mean coverage.
- bowtie2coverage.py: python script to convert bowtie1 output to per position
  coverage table. NOT used in the final analysis
- cap3\_all.job: Slurm script to do cap3 assembly with the default parameters
- cat\_dna\_rna.sh: bash script to combine edited dna and rna sequences
- de-edit\_transcripts.py: alternative way to de-edit transcript sequences
- dist\_calc.py: calculate p-distances for an alignment and print basic
  statistics
- do4all: bash script to run the same command on a set of files. Option to save
  stout to a file.
- etandem2gff3.py: convert etandem output to a gff file. Options to output gff
  for repeat clusters or individual repeat units.
- etandem2gff3\_old.py: OLD version of the script
- extract\_by\_length.job: Slurm script to filter sequences by length. NOT used
  in the final analysis 
- extract\_features.py: extract sequences based on blast hits
- extract\_repeat\_chunks.py: extract repeat units of a certain size range based on gff file 
- extract\_repeat\_chunks\_bak.py: extract repeat clusters rather than individual repeats. NOT usef in the final analysis
- extract\_repeat\_chunks\_old.py: OLD version of the script
- extract\_transcripts.sh: bash script to extract gene-specific transcripts
  based on blast outfile
- extract\_transcripts\_old.sh: OLD version of the script
- gff2pdf.py: a file to create a pdf file based on gff file
- irf2gff.py: a script to convert irf output to a gff file. NOT used in the
  final analysis
- map\_subsequence\_to\_sequence.py: script to map subsequences (e.g.,
  transcripts) on a full sequence (e.g., chromosomes).
- polyAT\_to\_gff.py: detect polyA/T runs and output their positions as a gff
  file
- process\_sra.job: SLURM script to submit process\_sra.sh script
- process\_sra.sh: bash script to process sra files
- run\_etandem.sh: bash script to run etandem command on multiple sequences
- tRNAscan.job: SLURM script to submit tRNAscan.sh script
- tRNAscan.sh: bash script to run tRNAscan command
- update\_gff\_etandem\_strand\_from\_blast.py: add repeat direction based on
  blast results to individual repeats. NOT used in the final analysis.
- update\_gff\_etandem\_strand\_from\_cluster\_blast.py: add repeat direction
  based on blast results to repeat clusters.
