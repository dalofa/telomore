# TELOMORE
Telomore is a tool for identifying and extracting telomeric sequences from **Oxford Nanopore** or **Illumina** sequencing reads of *Streptomyces* that have been excluded from a de novo assembly. It processes sequencing data to extend assemblies, generate quality control (QC) reports, and produce finalized assemblies with the telomere included.

## Before running Telomore
Telomore does not identify linear contigs but rather rely on the user to provide that information in
the header of the fasta-reference file. 

## Usage
```
python [telomore.py] --mode <mode> --reference <reference.fasta> [options]
```

Required Arguments
--mode: Specify the sequencing platform. Options: nanopore or illumina.
--reference: Path to the reference genome file in FASTA format.

Nanopore-Specific Arguments
--single: Path to a single gzipped FASTQ file containing Nanopore reads.

Illumina-Specific Arguments
--read1: Path to gzipped FASTQ file for Illumina read 1.
--read2: Path to gzipped FASTQ file for Illumina read 2.

Optional Arguments
--threads: Number of threads to use (default: 1).
--keep: Retain intermediate files (default: False).
--quiet: Suppress console logging.

## Process overview
The process is as follows:
1. **Map Reads:**
Reads are mapped against all contigs in a reference using either minimap2 or Bowtie2.
2. **Extract Extending Reads**
Extending reads that are mapped to the ends of linear contigs are extracted.
3. **Build Consensus**
The terminal extending reads from each end is used to construct a consensus using either lamassemble or mafft + EMBOSS cons
4. **Align and Attach consensus**
The consensus for each end is aligned to the reference and used to extend it.
5. **Trim Extended Replicon**
In a final step, all terminally mapped reads are mapped to the new extended reference and used to trim away spurious sequence, based on read-support.

Output
Extended Assembly: Written to <basename>.02.trimmed.fasta.
QC Reports: Saved in a folder named <basename>_seqtype_QC.
Logs: Written to telomore.log and <basename>.seqtype.cons.log.txt.


## Dependencies (CLI-tools)
* Minimap, version 2.25 or higher
* Bowtie2
* Samtools
* Lamassemble
* LAST-DB
* Lamassemble
* Mafft
* Emboss tools (cons specifically)

These can be installed using conda
```
conda create -n telo_env \
minimap2=2.25 \
bowtie2 \
samtools \
lamassemble \
last \
mafft \
emboss \
-y
```





### Python Dependencies
Python packages nessesary to run the script (generated using pipreqs):
* Bio==1.8.0
* biopython==1.78
* GitPython==3.1.44
* pysam==0.23.3

