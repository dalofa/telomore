"""Functions for generating useful QC metrics from the telomore script."""

import csv
import os
import tempfile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam

from .cmd_tools import map_and_sort, map_and_sort_illumina
from .fasta_tools import (
    cat_and_derep_fastq,
    check_fastq_order,
    dereplicate_fastq,
    merge_fasta,
)
from .map_tools import sam_to_fastq, sam_to_readpair


def qc_map(
    extended_assembly: str, left: str, right: str, output_handle: str, t: int = 1
) -> None:
    """
    Generate QC alignment of terminal reads against extended assembly (Nanopore).

    Collects terminal reads from left and right SAM files, converts them to
    FASTQ, deduplicates, and maps them back to the extended assembly. Used
    to validate the quality of the extension by visualizing read support.

    Parameters
    ----------
    extended_assembly : str
        Path to FASTA file of extended assembly with consensus attached
    left : str
        Path to SAM file containing left-terminal reads
    right : str
        Path to SAM file containing right-terminal reads
    output_handle : str
        Path for output sorted BAM file with QC alignments
    t : int, default=1
        Number of threads for mapping

    Returns
    -------
    None
        Writes QC alignment BAM to output_handle

    Notes
    -----
    QC mapping workflow:
    1. Creates temporary FASTQ file
    2. Converts left SAM reads to FASTQ (appends to temp file)
    3. Converts right SAM reads to FASTQ (appends to temp file)
    4. Deduplicates the combined FASTQ to remove redundant reads
    5. Maps deduplicated reads to extended assembly using minimap2
    6. Sorts and indexes the resulting BAM file
    7. Cleans up temporary FASTQ file

    The resulting BAM file can be visualized in IGV or analyzed with
    samtools to assess:
    - Coverage across extended regions
    - Read support for consensus sequences
    - Consistency of read alignments at telomeres

    Deduplication is critical because the same read may appear in both
    left and right terminal sets if it maps near both ends.
    """
    # The file has to be mode=w to create the correct type of object
    # for sam_to_fastq
    with tempfile.NamedTemporaryFile(
        suffix='.fastq', delete=False, mode='w'
    ) as temp_fastq:
        temp_fastq_path = temp_fastq.name
        sam_to_fastq(left, temp_fastq)
        sam_to_fastq(right, temp_fastq)
    dereplicate_fastq(fastq_in=temp_fastq_path, fastq_out=temp_fastq_path)
    map_and_sort(extended_assembly, temp_fastq_path, output_handle, t)
    os.remove(temp_fastq_path)


def qc_map_illumina(
    extended_assembly: str,
    left_sam: str,
    right_sam: str,
    fastq_in1: str,
    fastq_in2: str,
    output_handle: str,
    t: int = 1,
) -> None:
    """
    Generate QC alignment of terminal reads against extended assembly (Illumina).

    Collects complete paired-end reads for all terminal alignments from left
    and right SAM files, deduplicates, and maps them back to the extended
    assembly. Preserves read pairing for accurate Illumina QC assessment.

    Parameters
    ----------
    extended_assembly : str
        Path to FASTA file of extended assembly with consensus attached
    left_sam : str
        Path to SAM file containing left-terminal read alignments
    right_sam : str
        Path to SAM file containing right-terminal read alignments
    fastq_in1 : str
        Path to original R1 FASTQ file (gzip compressed)
    fastq_in2 : str
        Path to original R2 FASTQ file (gzip compressed)
    output_handle : str
        Path for output sorted BAM file with QC alignments
    t : int, default=1
        Number of threads for mapping

    Returns
    -------
    None
        Writes QC alignment BAM to output_handle

    Raises
    ------
    Exception
        If FASTQ files are not properly paired or ordered

    Notes
    -----
    QC mapping workflow for paired-end data:
    1. Creates temporary directory for intermediate files
    2. Extracts both R1 and R2 for reads in left SAM from original FASTQs
    3. Extracts both R1 and R2 for reads in right SAM from original FASTQs
    4. Concatenates and deduplicates R1 files separately
    5. Concatenates and deduplicates R2 files separately
    6. Validates that R1 and R2 files are properly synchronized
    7. Maps paired reads to extended assembly using BWA-MEM
    8. Sorts and indexes the resulting BAM file
    9. Cleans up temporary directory

    The paired-end approach ensures:
    - Proper insert size analysis for extended regions
    - Better mapping quality through paired information
    - Accurate assessment of consensus support from both read ends

    File order validation is critical - BWA-MEM requires synchronized
    R1/R2 pairs, and the check prevents mapping with mismatched pairs.
    """
    # get left paired read
    with (
        tempfile.TemporaryDirectory()
    ) as temp_dir:  # ensures files are deleted after usage
        # Create multiple temporary files in the temporary directory
        l_tmp1 = os.path.join(temp_dir, 'terminal_left_reads_1.fastq')
        l_tmp2 = os.path.join(temp_dir, 'terminal_left_reads_2.fastq')
        r_tmp1 = os.path.join(temp_dir, 'terminal_right_reads_1.fastq')
        r_tmp2 = os.path.join(temp_dir, 'terminal_right_reads_2.fastq')
        a_tmp1 = os.path.join(temp_dir, 'all_terminal_reads_1.fastq')
        a_tmp2 = os.path.join(temp_dir, 'all_terminal_reads_2.fastq')

        sam_to_readpair(
            sam_in=left_sam,
            fastq_in1=fastq_in1,
            fastq_in2=fastq_in2,
            fastq_out1=l_tmp1,
            fastq_out2=l_tmp2,
        )

        # get right paired read

        sam_to_readpair(
            sam_in=right_sam,
            fastq_in1=fastq_in1,
            fastq_in2=fastq_in2,
            fastq_out1=r_tmp1,
            fastq_out2=r_tmp2,
        )

        # collect the paired read files:
        cat_and_derep_fastq(fastq_in1=l_tmp1, fastq_in2=r_tmp1, fastq_out=a_tmp1)

        cat_and_derep_fastq(fastq_in1=l_tmp2, fastq_in2=r_tmp2, fastq_out=a_tmp2)

        if check_fastq_order(a_tmp1, a_tmp2):
            map_and_sort_illumina(
                reference=extended_assembly,
                read1=a_tmp1,
                read2=a_tmp2,
                output=output_handle,
                threads=t,
            )
        else:
            raise Exception('FASTQ files are not properly paired or ordered.')


def cons_genome_map(
    left_cons: str,
    right_cons: str,
    polished_genome: str,
    output_handle: str,
    t: int = 1,
) -> None:
    """
    Map consensus sequences against the polished reference genome.

    Merges left and right consensus sequences and aligns them to the final
    polished genome to verify their placement and identify any issues with
    consensus quality or positioning. Used for QC validation of consensus.

    Parameters
    ----------
    left_cons : str
        Path to FASTA file containing left consensus sequence
    right_cons : str
        Path to FASTA file containing right consensus sequence
    polished_genome : str
        Path to FASTA file of polished/final genome
    output_handle : str
        Path for output sorted BAM file with consensus alignments
    t : int, default=1
        Number of threads for mapping

    Returns
    -------
    None
        Writes consensus alignment BAM to output_handle

    Notes
    -----
    This QC mapping helps identify:
    - Whether consensus sequences map uniquely to their expected locations
    - If consensus contains repeats that map to multiple locations
    - Quality of consensus alignment (mismatches, soft-clipping)
    - Whether consensus extends correctly from the reference

    The temporary merged FASTA file 'all_cons.fasta' is created in the
    current directory and not automatically cleaned up.

    Consensus sequences should map with high identity to their respective
    ends. Multiple mappings or poor alignment quality suggests the consensus
    may not represent true telomeric extension.
    """
    merge_fasta(left_cons, right_cons, 'all_cons.fasta')
    map_and_sort(polished_genome, 'all_cons.fasta', output_handle, t)


def cons_cons_map(
    left_cons: str, right_cons: str, output_handle: str, t: int = 1
) -> None:
    """
    Map left consensus against right consensus to detect similarity.

    Aligns left and right consensus sequences against each other to identify
    potential circularization or repetitive telomeric sequences. If consensus
    sequences map to each other, it may indicate the chromosome is circular
    or contains telomeric repeats.

    Parameters
    ----------
    left_cons : str
        Path to FASTA file containing left consensus (used as reference)
    right_cons : str
        Path to FASTA file containing right consensus (used as query)
    output_handle : str
        Path for output sorted BAM file with cross-consensus alignments
    t : int, default=1
        Number of threads for mapping

    Returns
    -------
    None
        Writes cross-consensus alignment BAM to output_handle

    Notes
    -----
    Interpretation of results:
    - No alignment: Linear chromosome with distinct telomeres (expected)
    - High-quality alignment: May indicate:
      * Circular chromosome where ends should connect
      * Telomeric repeat arrays present at both ends
      * Potential artifact if sequences shouldn't match

    This QC check is particularly useful for:
    - Bacterial genomes where circularity is expected
    - Identifying repetitive telomeric sequences
    - Validating that linear chromosome ends are truly distinct

    Maps right consensus (query) against left consensus (reference) using
    minimap2 single-read mode with sorting and indexing.
    """
    map_and_sort(left_cons, right_cons, output_handle, t)


def cons_length(cons_file: str, output_handle: str, offset: int = 100) -> None:
    """
    Write consensus sequence length statistics to TSV file.

    Calculates and records the length of consensus sequences with and without
    an offset adjustment. The offset represents the amount of original reference
    sequence included in the consensus file for context.

    Parameters
    ----------
    cons_file : str
        Path to FASTA file containing consensus sequences
    output_handle : str
        Path for output TSV file with length statistics
    offset : int, default=100
        Number of bases of original reference included in consensus

    Returns
    -------
    None
        Writes TSV with columns: seq_id, end_cons, full_cons

    Notes
    -----
    Output TSV format:
    - Header: seq_id, end_cons, full_cons
    - seq_id: Identifier of the consensus sequence
    - end_cons: Length of true extension (full_cons - offset)
    - full_cons: Total length including offset region

    The offset adjustment is important because consensus building may
    include some bases from the original reference sequence for context
    and alignment purposes. The 'end_cons' value represents only the
    novel sequence extending beyond the original assembly.

    Example:
    - Full consensus: 250bp
    - Offset: 100bp
    - True extension: 150bp (reported as end_cons)

    Used for summarizing extension results across multiple contigs.
    """
    cons_file = SeqIO.parse(cons_file, 'fasta')
    header = ['seq_id', 'end_cons', 'full_cons']
    tsv_log = []
    tsv_log.append(header)

    for record in cons_file:
        seq_id = record.id
        seq_len = len(record)
        gen_len = int(seq_len) - offset
        tsv_log.append([seq_id, gen_len, seq_len])

    with open(output_handle, 'w', newline='') as tsv_file:
        writer = csv.writer(tsv_file, delimiter='\t')
        writer.writerows(tsv_log)


def map_to_depth(bam_file: str, output_handle: str) -> None:
    """
    Generate position-by-position depth of coverage from BAM file.

    Extracts coverage depth at every position in the reference using samtools
    depth. Creates a tab-delimited file showing reference name, position, and
    coverage at each position. Used for visualizing coverage profiles.

    Parameters
    ----------
    bam_file : str
        Path to input BAM alignment file
    output_handle : str
        Path for output depth file

    Returns
    -------
    None
        Writes depth information to output_handle

    Notes
    -----
    Uses 'samtools depth -aa' which:
    - -aa: Output absolutely all positions, including zero-coverage
    - Ensures complete coverage profile even for uncovered regions

    Output format (tab-delimited):
    - Column 1: Reference sequence name
    - Column 2: Position (1-based)
    - Column 3: Coverage depth at that position

    The output file can be:
    - Plotted to visualize coverage across genome
    - Used to identify low-coverage regions
    - Analyzed to assess quality of consensus extensions
    - Imported into visualization tools like R or Python

    Particularly useful for QC visualization to show how coverage
    changes across telomeric regions and consensus extensions.
    """
    pysam.depth('-aa', bam_file, '-o', output_handle)


def finalize_log(log: str, right_fasta: str, left_fasta: str) -> None:
    """
    Finalize extension log by prepending final consensus lengths and sequences.

    Rewrites the log file with a summary section at the top showing the final
    validated consensus lengths after trimming, followed by the original log
    content documenting the extension process. Extracts trimmed consensus
    sequences and includes them in the final summary.

    Parameters
    ----------
    log : str
        Path to extension log file to finalize (will be overwritten)
    right_fasta : str
        Path to FASTA file containing right consensus sequence
    left_fasta : str
        Path to FASTA file containing left consensus sequence

    Returns
    -------
    None
        Overwrites log file with finalized version including summary header

    Notes
    -----
    Log processing steps:
    1. Reads existing log content
    2. Extracts original consensus lengths from line 4
    3. Extracts trimming information from last two lines
    4. Calculates final lengths: original_length - trimmed_bases
    5. Extracts trimmed portion of consensus sequences
    6. Writes new header section with final results
    7. Appends original log content below header

    Final log structure:
    - FINAL GENOME EXTENSION header
    - Final consensus lengths (may show 'rejected' if validation failed)
    - Final consensus sequences (trimmed portions only)
    - Separator line
    - Original log content (initial consensus, trimming details)
    - Closing separator

    Handles rejected consensus:
    - If 'rejected' in trim log: Shows 'rejected' instead of length
    - Creates empty SeqRecord for rejected consensus
    - For accepted consensus: Shows length and trimmed sequence

    The final log provides a complete record of:
    - What extensions were added (top section)
    - How they were generated and validated (original log below)
    """
    file = open(log)
    log_cont = file.readlines()
    file.close()
    # Org lengths of consensus added
    length_lines = log_cont[3]
    left_len = int(length_lines.split('\t')[0].split(':')[1])
    right_len = int(length_lines.split('\t')[1].split(':')[1])

    # get the number of bases trimmed off
    trim_left = log_cont[-2].split(' ')[-1]
    trim_right = log_cont[-1].split(' ')[-1]
    left_seq = SeqIO.read(left_fasta, 'fasta')
    right_seq = SeqIO.read(right_fasta, 'fasta')

    if trim_left.rstrip() == 'rejected':
        new_left = 'rejected'
        left_seq = SeqRecord(Seq(''))
    else:
        new_left = left_len - int(trim_left)
        left_seq = left_seq[int(trim_left) :]
    if trim_right.rstrip() == 'rejected':
        new_right = 'rejected'
        right_seq = SeqRecord(Seq(''))
    else:
        new_right = right_len - int(trim_right)
        right_seq = right_seq[0:new_right]

    final_lengths = 'left_cons:{}\tright_consensus:{}'.format(new_left, new_right)

    # write to log file
    file = open(log, 'w')
    file.write(
        '=============================================================================='
    )
    file.write('\nFINAL GENOME EXTENSION')
    file.write(
        '\n==============================================================================\n'
    )
    file.write(final_lengths)
    file.write('\n>left_cons\n')
    file.write(str(left_seq.seq))
    file.write('\n>right_cons\n')
    file.write(str(right_seq.seq))
    file.write('\n')
    for line in log_cont:
        file.write(line)
    file.write(
        '==============================================================================\n'
    )
    file.close()
