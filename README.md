# TELOMORE
A pipeline to reconstitute telomeric sequences in *Streptomyces* genomes excluded during de novo assembly from either Oxford Nanopore or Illumina sequencing data. In brief, the tool looks for reads that extend an assembly, builds a consensus of the extending sequence and attaches it to the assembly.

## Usage
```
usage: telomore.py [-h] -f FASTQ -r REFERENCE -m {nanopore,illumina} [-t THREADS] [-k]

Recover potential telomeric sequences from Streptomyces genomes using Oxford Nanopore or illumina sequence data.
- OUTPUT: An extended assembly is written to basename.02.trimmed.fasta and QC-maps are written to a folder
  named basename_seqtype_QC.
- LOG: A run-log is written to telomore.log and a result-log is written to basename.seqtype.cons.log.txt.

options:
  -h, --help            show this help message and exit
  -f FASTQ, --fastq FASTQ
                        Path to gzipped fastq-file
  -r REFERENCE, --reference REFERENCE
                        Path to reference file (.fasta, .fna, or .fa)
  -m {nanopore,illumina}, --mode {nanopore,illumina}
                        Choose which mode to run
  -t THREADS, --threads THREADS
                        Threads to use. Default is 1
  -k, --keep            Flag to keep intermediate files. Default is False
```


## Process overview
The process is as follows:

0. **Identify the chromosome:**
The longest contig is extracted and assumed to be the chromosome.
1. **Map reads to chromosome**
Nanopore reads are mapped to the chromosome using minimap2, while Illumina reads are mapped using Bowtie2. Terminal extending reads are extracted from this maps.
2. **Generate consensus:**
A consensus are constructed for each terminus of the genome using the extending reads. This is done using lamassemble for nanopore reads or mafft for illumina reads.
3. **Extend assembly with consensus:**
Each terminus of the chromsome is extended with corresponding consensus by mapping this consensus onto the chromsome.
4. **Trimming of newly attached sequence:**
The previously identified terminal reads are mapped to extended assembly and the sequence is trimmed based on read-support as specified in the result log.
5. **QC and clean-up:**
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

