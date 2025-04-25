# TELOMORE
A pipeline to reconstitute telomeric sequences in *Streptomyces* genomes excluded during de novo assembly from either Oxford Nanopore or Illumina sequencing data. In brief, the tool looks for reads that extend an assembly, builds a consensus of the extending sequence and attaches it to the assembly.

Telomore does not check contigs for circularity but instead relies on the user to tag contigs as "linear" in the fasta_header.

## Usage
```
usage: telomore.py [-h] -m {nanopore,illumina} [--single SINGLE] [--read1 READ1] [--read2 READ2] -r REFERENCE [-t THREADS] [-k]

Recover potential telomeric sequences from Streptomyces genomes using Oxford Nanopore or illumina sequence data.

- INPUT: Takes a fasta file with linear contigs tagged "linear" in the fasta_header and sequencing reads in fq.gz format.
- OUTPUT: An extended assembly is written to basename.02.trimmed.fasta and QC-maps are written to a folder
  named basename_seqtype_QC.
- LOG: A run-log is written to telomore.log and a result-log is written to basename.seqtype.cons.log.txt.

options:
  -h, --help            show this help message and exit
  -m {nanopore,illumina}, --mode {nanopore,illumina}
                        Choose which mode to run.
                                --mode=nanopore takes a single read-file, specified using --single
                                --mode=illumina-mode takes two read-files, specified using --read1 and --read2
  --single SINGLE       Path to a single gzipped nanopore fastq-file
  --read1 READ1         Path to gzipped illumina read1 fastq-file
  --read2 READ2         Path to gzipped illumina read2 fastq-file
  -r REFERENCE, --reference REFERENCE
                        Path to reference file (.fasta, .fna, or .fa)
  -t THREADS, --threads THREADS
                        Threads to use. Default is 1
  -k, --keep            Flag to keep intermediate files. Default is False

```


## Process overview
The process is as follows:
1. **Map Reads:**
Reads are mapped against all contigs in a reference using either minimap2 or Bowtie2.
2. **Extract Extending Reads**
Extending reads that are mapped to the ends of linear contigs are extracted.
3. **Build Consensus**
The terminal extending reads from each end is used to construct a consensus using either lamassemble [] or mafft [] and EMBOSS cons [].
4. **Align and Attach consensus**
The consensus for each end is aligned to the reference and used to extend it.
5. **Trim Extended Replicon**
In a final step, all terminally mapped reads are mapped to the new extended reference and used to trim away spurious sequence, based on read-support.

Output:
A number of QC-maps are generated and moved to a folder called *_seqtype_QC, where * is the basename of the reference supplied.
In here three bam files can be found together with their references:
- 1: The terminal reads mapped against the genome+full consensus. These are denoted with *.nontrimmed*
- 2: The terminal reads mapped against the final genome+trimmed_consensus. These are denoted with *.trimmed.*
- 3: The two consensuses mapped against the genome



## Dependencies (CLI-tools)
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

