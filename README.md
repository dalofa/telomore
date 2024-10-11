# TELOMORE
A pipeline for fishing out telomores trimmed off a linear genome (Streptomyces) by Flye.

## Process overview


![telo_crop_pipeline](https://github.com/dalofa/telomore/assets/83669966/1cdc3192-3821-4b70-8767-8613accaced0)

The process is as follows:

0. Identify the genome:
The longest contig is extracted and assumed to be the chromosome.
1. Mapping of reads and extraction of soft-clipped reads:
Nanopore reads are mapped to the genome using minimap2 and terminally mapped reads that are soft-clipped are extracted.
3. Consensus building and attachment:
Using LAST and lamassemble a consensus is constructed from soft-clipped sequence that extends the genome + 500 bases. The resulting consensus
is aligned to teh chromosome and any sequence that extend beyond the chromosome ends are added, creating a temporary consensus+genome+consensus.
4. Trimming of newly attached sequence:
The previously identified terminal reads are mapped to genome+consensus and the sequence is trimmed. This is done by
removing bases from the ends until at least 5 reads support the position and at least 70% of the reads agree with the position.
5. QC and clean-up:
A number of QC-maps are generated and moved to a folder called *_QC, where * is the basename of the reference supplied.
In here three bam files can be found together with their references:
- 1: The terminal reads mapped against the genome+full consensus. These are denoted with *.nontrimmed*
- 2: The terminal reads mapped against the final genome+trimmed_consensus. These are denoted with *.trimmed.*
- 3: The two consensuses mapped against the genome

#### Matching QC mapping files to fasta files:

A number prefix is used to match a reference to a given bam file, such that BASE.02.reads.trimmed.bam is a mapping of reads
onto the BASE.02.trimmed.cons.fasta file. In the same vein BASE.02.cons.trimmed.bam is a mapping of the consensus seq onto
the BASE.02.trimmed.cons.fasta file.

A file called *_consensus.log.txt can be found supplying information about the trimming of the attached sequence.

The final file of genome+trimmed consensus is called *.trimmed.cons.fasta


#### BE AWARE TUE!
Currently the script also leaves the map used to trim the consensus sequence + the matching fasta file, which were the maps we used to QC the trimming.
The map is named *.nontrimmed.map.sam.sort.bam and the corresponding fasta file is *_genome.stitch.fasta (which contains the genome + untrimmed consensus)


## Running on NBC-shared
Before running the script lamassemble and last must be loaded.
> module load lamassemble/1.4.2 last/1448 minimap2/2.25

To run the script:
>python3 path_to_telomore/main.py -f fastq -r reference -t threads

More details:
```
usage: main.py [-h] -f FASTQ -r REFERENCE [-t THREADS] [-k]

Recover potential telomeric seq from Streptomyces Oxford Nanopore data

options:
  -h, --help            show this help message and exit
  -f FASTQ, --fastq FASTQ
                        Path to gzipped fastq-file
  -r REFERENCE, --reference REFERENCE
                        Path to gzipped reference fastq-file
  -t THREADS, --threads THREADS
                        Threads to use. Default is 1
  -k, --keep            Flag to keep intermediate files. Default is False
```
### Dependencies
Minimap, version 2.25 or higher
Bowtie2
Samtools
Lamassemble
LAST-DB
Lamassemble
Mafft
Emboss tools (cons specifically)



### Python Dependencies
Python packages nessesary to run the script (generated using pipreqs):
Bio==1.7.0
biopython==1.81
pysam==0.21.0

