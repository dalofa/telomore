# Telomore: Project Overview

## Purpose

Telomore is a bioinformatics tool designed to recover and extend telomeric sequences from Streptomycetes genome assemblies. When performing *de novo* genome assembly, telomeric sequences at chromosome ends are often excluded or incomplete due to low coverage or repetitive structures. Telomore addresses this limitation by:

1. Identifying reads that extend beyond the ends of linear contigs in an assembly
2. Building consensus sequences from these extending reads
3. Attaching and validating these consensus sequences to extend the assembly
4. Trimming spurious bases based on read support

The tool supports both **Oxford Nanopore** long-read sequencing and **Illumina** short-read paired-end sequencing platforms, with platform-specific algorithms optimized for each technology's characteristics.

## Project Structure

```text
telomore/
├── src/telomore/
│   ├── __init__.py              # Package initialization
│   ├── app.py                   # Main application entry point and workflow orchestration
│   └── utils/
│       ├── __init__.py          # Utils package initialization
│       ├── arg_parser.py        # Command-line argument parsing and logging setup
│       ├── cmd_tools.py         # External tool wrappers (minimap2, bowtie2, samtools, etc.)
│       ├── fasta_tools.py       # FASTA/FASTQ file manipulation utilities
│       ├── map_tools.py         # Read mapping analysis and consensus extension logic
│       ├── qc_reports.py        # Quality control report generation
│       └── classes_and_small_func.py  # Replicon class for managing file paths
├── pyproject.toml               # Project configuration and dependencies
├── environment.yml              # Conda environment specification
├── README.md                    # User documentation
└── CITATION.cff                 # Citation metadata

```

## Architecture and Implementation

### Core Workflow (app.py)

The `main()` function in `app.py` orchestrates the entire telomere extension workflow:

1. **Dependency checking**: Verifies all required external tools are installed
2. **Input validation**: Identifies linear contigs tagged with "linear" in FASTA headers
3. **Read mapping**: Maps all reads to the reference assembly
4. **Terminal read extraction**: Identifies reads extending beyond contig ends
5. **Consensus generation**: Builds consensus sequences from extending reads
6. **Assembly extension**: Aligns and attaches consensus to contig ends
7. **Consensus trimming**: Validates and trims consensus based on read support
8. **QC generation**: Creates quality control BAM files for manual inspection
9. **Output finalization**: Combines extended contigs with unchanged contigs

### Key Components

#### 1. Replicon Class (classes_and_small_func.py)

The `Replicon` class encapsulates all file paths associated with processing a single linear contig:

- Input files: Original FASTA, mapped reads
- Filtered reads: Left and right terminal reads in SAM/FASTQ formats
- Consensus files: Left and right consensus sequences and alignments
- Extension files: Truncated references, consensus maps, extended assemblies
- Output files: Trimmed assemblies, QC BAM files, logs

This design provides clean separation of concerns and makes file management straightforward.

#### 2. Command-Line Tools Integration (cmd_tools.py)

Telomore wraps several external bioinformatics tools:

**For Nanopore reads:**

- `minimap2`: Long-read alignment with parameters optimized for soft-clipping
- `lamassemble`: Consensus generation for dissimilar reads
- `lastdb`/`last-train`: Training alignment parameters on the reference

**For Illumina reads:**

- `bowtie2`: Short-read paired-end alignment with local mode
- `mafft`: Multiple sequence alignment for similar reads
- `cons` (EMBOSS): Consensus calling from aligned sequences

**Shared tools:**

- `samtools`: SAM/BAM manipulation, sorting, and indexing

Each wrapper function includes:

- Input validation (file existence, parameter types)
- Comprehensive logging at each step
- Error handling with detailed error messages
- Automatic cleanup of temporary files
- Thread-safe parallel execution support

#### 3. FASTA/FASTQ Utilities (fasta_tools.py)

Provides utilities for sequence file manipulation:

- `get_linear_elements()`: Parses FASTA headers to identify linear contigs
- `extract_contig()`: Extracts individual contigs from multi-FASTA files
- `strip_fasta()`: Removes bases from sequence ends to avoid alternative mapping sites
- `build_extended_fasta()`: Reconstructs the final assembly with extended contigs
- `dereplicate_fastq()`: Removes duplicate reads by ID
- `check_fastq_order()`: Validates paired-end read synchronization

#### 4. Mapping Analysis (map_tools.py)

Contains the core logic for telomere extension:

**Terminal read identification:**

- `get_terminal_reads()`: Extracts reads mapping to first/last 20bp of contigs
- `get_left_soft()`/`get_right_soft()`: Filters soft-clipped reads extending beyond genome edges
- Uses regex parsing of CIGAR strings to identify soft-clipped regions
- Deduplicates by keeping the alignment with most mapped bases

**Consensus extension logic (`stitch_telo()`):**

- Analyzes soft-clipped portions of consensus alignments to reference
- Extracts only the overhanging (extending) portion of consensus sequences
- Handles edge cases: empty consensus, unmapped consensus, non-extending consensus
- Concatenates: left_consensus + original_contig + right_consensus

**Consensus validation and trimming:**

- `trim_by_map()`/`trim_by_map_illumina()`: Trims consensus based on read support
- `get_support_info()`: Calculates per-position coverage and base matching
- Trimming criteria: coverage ≥ threshold AND match ratio ≥ 0.7
- Default thresholds: Nanopore (cov≥5, Q≥10), Illumina (cov≥1, Q≥30)

**Key algorithms:**

- `mapped_bases()`: Parses CIGAR strings to count reference-consuming operations (M, D, N, X, =)
- `cigar_maps_more_bases()`: Compares alignments to select best mapping per read
- `revcomp_reads()`: Reverse-complements reads and quality scores for left-side processing

#### 5. Quality Control (qc_reports.py)

Generates validation files for manual inspection:

- `qc_map()`: Maps terminal reads to extended assembly for Nanopore
- `qc_map_illumina()`: Maps terminal read pairs to extended assembly for Illumina
- `finalize_log()`: Creates detailed logs showing bases added/trimmed with sequences
- Uses temporary files to collect and deduplicate reads before final mapping

#### 6. Argument Parsing (arg_parser.py)

- `get_args()`: Defines CLI interface with argparse
- Validates mode-specific requirements (single FASTQ for Nanopore, paired FASTQ for Illumina)
- `setup_logging()`: Configures dual output to console and file
- Provides helpful error messages for missing required parameters

## Algorithm Details

### Why Streptomycetes-Specific?

Streptomycetes genomes have **Terminal Inverted Repeats (TIRs)** - identical sequences at both chromosome ends. To prevent consensus sequences from incorrectly mapping to the opposite end, Telomore:

1. Truncates half of each contig before consensus mapping
2. Maps left consensus only to left-truncated reference
3. Maps right consensus only to right-truncated reference

This ensures consensus sequences can only extend their intended end.

### Nanopore vs. Illumina Strategies

**Nanopore (long reads, higher error rate):**

- Uses `lamassemble` with trained LAST parameters for consensus of dissimilar sequences
- Default coverage threshold = 5x (higher to account for errors)
- Default quality threshold = Q10 (lower due to platform characteristics)
- Single-pass mapping with soft-clipping to capture extensions

**Illumina (short reads, lower error rate):**

- Uses `mafft + cons` assuming reads are more similar
- Default coverage threshold = 1x (lower due to higher accuracy)
- Default quality threshold = Q30 (higher, leveraging platform accuracy)
- Requires read pair extraction and synchronization for mapping

### Soft-Clipping Strategy

Soft-clipping occurs when a read extends beyond the reference. Example:

```text
Reference:  ATCGATCG|
Read:       ATCGATCGATCG
CIGAR:      8M4S (8 matched, 4 soft-clipped)
```

Telomore:

1. Identifies reads with soft-clips extending beyond reference boundaries
2. Extracts only the soft-clipped portions (the "extending" bases)
3. Builds consensus from these extending portions
4. Maps consensus back to validate it genuinely extends the assembly

### Consensus Trimming Algorithm

After generating an extended assembly, Telomore validates it by:

1. Mapping terminal reads to the extended sequence
2. Calculating per-position coverage and base match ratios
3. Iterating from the consensus ends inward until finding positions meeting criteria:
   - Coverage ≥ threshold
   - Match ratio ≥ 0.7 (70% of reads support the reference base)
   - Base quality ≥ threshold

This trimming removes spurious bases while keeping well-supported extensions.

## File Formats and Outputs

### Input Requirements

1. **Reference FASTA**: Contigs with "linear" keyword in description line

   ```text
   >contig_1 linear
   ATCGATCG...
   ```

2. **Sequencing reads**:
   - Nanopore: Single gzipped FASTQ file
   - Illumina: Two gzipped paired-end FASTQ files

### Output Structure

```text
{reference_basename}_{np|ill}_telomore/
├── {contig}_telomore_extended.fasta      # Final extended contig
├── {contig}_telomore_untrimmed.fasta     # Pre-trimming extended contig
├── {contig}_telomore_QC.bam              # Terminal reads mapped to final
├── {contig}_telomore_QC.bam.bai          # BAM index
├── {contig}_telomore_ext_{np|ill}.log    # Detailed extension log
└── {reference_basename}_telomore.fasta   # All contigs (extended + unchanged)
```

### Log File Format

The extension logs document:

1. **FINAL GENOME EXTENSION**: Actual bases added after trimming
2. **INITIAL CONSENSUS**: Raw consensus lengths before validation
3. **CONSENSUS TRIMMING**: Trimming rules and bases removed

## Dependencies and External Tools

Telomore orchestrates multiple specialized bioinformatics tools rather than reimplementing complex algorithms:

- **minimap2**: State-of-the-art long-read aligner
- **bowtie2**: Widely-used short-read aligner
- **samtools**: Standard for SAM/BAM manipulation
- **lamassemble**: Overlap-layout-consensus assembler for dissimilar sequences
- **mafft**: Multiple sequence alignment
- **EMBOSS cons**: Consensus calling
- **LAST**: Alignment with trainable parameters

Python dependencies:

- **BioPython**: Sequence parsing and manipulation
- **pysam**: Python interface to samtools

## Design Principles

1. **Modularity**: Clear separation between CLI tools, file I/O, and analysis logic
2. **Reproducibility**: Comprehensive logging of all steps and parameters
3. **Validation**: QC files enable manual verification of extensions
4. **Flexibility**: User-configurable thresholds for different sequencing depths
5. **Safety**: Only extends contigs explicitly tagged as "linear"
6. **Cleanup**: Automatic removal of intermediate files (unless `--keep` specified)

## Limitations and Considerations

- Requires pre-identification of linear contigs by the user
- Optimized for Streptomycetes (but applicable to other organisms)
- Depends on availability of reads extending beyond assembly
- Low coverage at chromosome ends may prevent extension
- Repetitive telomeric sequences may complicate consensus generation
