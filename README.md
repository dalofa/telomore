# TELOMORE
Telomore is a tool for identifying and extracting telomeric sequences from **Oxford Nanopore** or **Illumina** sequencing reads of Streptomycetes that have been excluded from a de novo assembly. It processes sequencing data to extend assemblies, generate quality control (QC) maps, and produce finalized assemblies with the telomere/recessed bases included.

## Before running Telomore
Telomore does not identify linear contigs but rather rely on the user to provide that information in
the header of the fasta-reference file. 

## Usage
```
telomore --mode <mode> --reference <reference.fasta> [options]
```

Required Arguments
- `--mode` Specify the sequencing platform. Options: nanopore or illumina.
- `--reference` Path to the reference genome file in FASTA format.

Nanopore-Specific Arguments
- `--single` Path to a single gzipped FASTQ file containing Nanopore reads.

Illumina-Specific Arguments
- `--read1` Path to gzipped FASTQ file for Illumina read 1.
- `--read2` Path to gzipped FASTQ file for Illumina read 2.

Optional Arguments
- `--coverage_threshold` Set the threshold for coverage to stop trimming during consensus trimming (Default is coverage=5 for ONT reads and coverage=1 for Illumina reads).
- `--quality_threshold` Set the Q-score required to count a read position in the coverage calculation during consensus trimming (Default is Q-score=10 for ONT reads and Q-score=30 for Illumina reads).
- `--threads` Number of threads to use (default: 1).
- `--keep` Retain intermediate files (default: False).
- `--quiet` Suppress console logging.

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

## Outputs
At the end of a run Telomore produces the following outputs:

```Output
├── {fasta_basename}_{seqtype}_telomore
│   ├── {contig_name}_telomore_extended.fasta
│   ├── {contig_name}_telomore_ext_{seqtype}.log
│   ├── {contig_name}_telomore_QC.bam
│   ├── {contig_name}_telomore_QC.bam.bai
│   ├── {contig_name}_telomore_untrimmed.fasta
│   └── {fasta_basename}_telomore.fasta
└── telomore.log # log containing run information.
```
In the folder there is a number of files generated for each contig considered:
| File Name | Description |
|-----------|-------------|
| `{contig_name}_telomore_extended.fasta` | Original contig sequence + added terminal bases - trimmed bases |
| `{contig_name}_telomore_ext_{seqtype}.log` | Log contianing information about bases added, trimmed off and final result. |
| `{contig_name}_telomore_QC.bam` | BAM file containing terminal reads mapped to `{contig_name}_telomore_extended.fasta`. Useful for manual inspection of the extension|
| `{contig_name}_telomore_QC.bam.bai` | Index file for the corresponding BAM file. |
| `{contig_name}_telomore_untrimmed.fasta` | Original contig sequence + added terminal bases |

Additionally, there is a fasta-file collecting all tagged linear contigs as they appear in `{contig_name}_telomore_extended.fasta` together with all non-linear contigs in the order they appear in the original file.

Inspecting the {contig_name}_QC.bam-file in IGV (Integrative Genomics Viewer) can be informative in evaluating the extended contig.

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
python=3.11.9 \
minimap2=2.25 \
bowtie2 \
samtools \
lamassemble \
last \
mafft \
emboss \
-y
```

This repo can then be downloaded using git clone, the conda enviroment activated and the tool installed
```
git clone https://github.com/dalofa/telomore
conda activate telo_env
cd telomore
pip install .
```

### Python Dependencies
Python packages nessesary to run the script (generated using pipreqs):
* Bio==1.8.0
* biopython==1.78
* GitPython==3.1.44
* pysam==0.23.3

