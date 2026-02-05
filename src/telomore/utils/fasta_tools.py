"""Utilities for handling fasta files."""

from itertools import zip_longest
import logging

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def check_fastq_order(file1: str, file2: str) -> bool:
    """
    Check if two FASTQ files have the same length and read order.

    Validates that paired-end FASTQ files are properly synchronized by ensuring
    they contain the same number of reads in the same order. This is critical
    for paired-end mapping tools which expect synchronized inputs.

    Parameters
    ----------
    file1 : str
        Path to first FASTQ file
    file2 : str
        Path to second FASTQ file

    Returns
    -------
    bool
        True if files are the same length with matching read IDs in order,
        False otherwise

    Notes
    -----
    This function:
    - Iterates through both files simultaneously using zip_longest
    - Compares read IDs at each position
    - Prints informative error message if mismatch found
    - Returns False immediately upon first mismatch

    A return value of False indicates the files cannot be used together
    for paired-end mapping without reordering or filtering.
    """
    handle1 = SeqIO.parse(file1, 'fastq')
    handle2 = SeqIO.parse(file2, 'fastq')

    # Iterate over reads, use zip_longest to not stop if one file is shorter than the other
    for i, (read1, read2) in enumerate(zip_longest(handle1, handle2)):
        if read1 is None or read2 is None:
            print(
                f'{file1} and {file2} are not the same length, diverging at read {i + 1}'
            )
            return False
        if read1.id != read2.id:
            print(
                f'Mismatch at read {i + 1} in files {file1} {file2}: {read1.id} != {read2.id}'
            )
            return False
    return True


def get_linear_elements(fasta_file: str) -> list[str]:
    """
    Extract contig names that are tagged as linear in a FASTA file.

    Parses a FASTA file to identify contigs with 'linear' in their description
    line. This is used to identify which contigs should be processed for
    telomere extension.

    Parameters
    ----------
    fasta_file : str
        Path to FASTA file where linear contigs are tagged

    Returns
    -------
    list of str
        List of contig IDs (record.id) for contigs with 'linear' in description

    Notes
    -----
    Expected FASTA header format for linear contigs:
        >contig_name linear
        or
        >contig_name [linear] some other description

    The 'linear' keyword can appear anywhere in the description line.
    Only the contig ID (before the first space) is returned, not the full
    description.

    Empty list is returned if no linear contigs are found, which causes
    the workflow to exit gracefully.
    """
    linear_list = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        if 'linear' in record.description:
            linear_list.append(record.id)
    return linear_list


def extract_contig(fasta_in: str, contig_name: str, fasta_out: str) -> None:
    """
    Extract a single contig from a multi-FASTA file.

    Searches through a FASTA file for a contig with the specified name and
    writes it to a new single-sequence FASTA file.

    Parameters
    ----------
    fasta_in : str
        Path to input multi-FASTA file
    contig_name : str
        Name of contig to extract (must match record.id exactly)
    fasta_out : str
        Path for output FASTA file containing only the extracted contig

    Returns
    -------
    None
        Writes extracted contig to fasta_out

    Notes
    -----
    - Only the first contig matching contig_name is extracted
    - If no match is found, no output file is created
    - The output FASTA retains the original sequence and description
    """
    for record in SeqIO.parse(fasta_in, 'fasta'):
        if record.id == contig_name:
            contig = record
            with open(fasta_out, 'w') as fq_file:
                SeqIO.write(sequences=contig, handle=fq_file, format='fasta')


def get_fasta_length(fasta_file: str, contig_name: str) -> int:
    """
    Get the sequence length of a specific contig in a FASTA file.

    Searches through a FASTA file for a contig with the specified name and
    returns its sequence length in bases.

    Parameters
    ----------
    fasta_file : str
        Path to FASTA file
    contig_name : str
        Name of contig whose length to retrieve (must match record.id exactly)

    Returns
    -------
    int
        Length of the contig sequence in bases

    Notes
    -----
    - Returns length of first matching contig
    - Returns None implicitly if contig not found (no explicit return statement)
    - Used to determine truncation boundaries for preventing alternative mappings
    """
    for record in SeqIO.parse(fasta_file, 'fasta'):
        if record.id == contig_name:
            length = len(record.seq)
            return length


def dereplicate_fastq(fastq_in: str, fastq_out: str) -> None:
    """
    Remove duplicate reads from a FASTQ file based on read ID.

    Creates a new FASTQ file containing only the first occurrence of each
    unique read ID. This prevents the same read from being counted multiple
    times in coverage calculations.

    Parameters
    ----------
    fastq_in : str
        Path to input FASTQ file (may contain duplicates)
    fastq_out : str
        Path for output deduplicated FASTQ file

    Returns
    -------
    None
        Writes deduplicated reads to fastq_out

    Notes
    -----
    - Deduplication is based solely on read.id (not sequence)
    - Order of first occurrences is preserved
    - Subsequent reads with the same ID are discarded
    - Useful when reads may map to multiple locations and appear in
      multiple SAM extractions
    """
    seen_reads = set()  # To store unique read identifiers and sequences
    unique_reads = []

    with open(fastq_in, 'r') as infile:
        for record in SeqIO.parse(infile, 'fastq'):
            # Dereplicate only on read_id to avoid identical reads with different mappings producing duplicates
            read_key = record.id

            if read_key not in seen_reads:
                seen_reads.add(read_key)
                unique_reads.append(record)

    with open(fastq_out, 'w') as outfile:
        SeqIO.write(unique_reads, outfile, 'fastq')


def cat_and_derep_fastq(fastq_in1: str, fastq_in2: str, fastq_out: str) -> None:
    """
    Concatenate two FASTQ files and remove duplicate reads.

    Combines two FASTQ files into a single output file and then removes
    duplicate reads based on read ID. This is useful for merging left and
    right terminal reads while ensuring each read appears only once.

    Parameters
    ----------
    fastq_in1 : str
        Path to first input FASTQ file
    fastq_in2 : str
        Path to second input FASTQ file
    fastq_out : str
        Path for output deduplicated FASTQ file

    Returns
    -------
    None
        Writes concatenated and deduplicated reads to fastq_out

    Notes
    -----
    This function operates in two stages:
    1. Concatenation: All reads from both input files are written to output
    2. Deduplication: The output file is overwritten with unique reads only

    The deduplication is performed by the dereplicate_fastq function, which
    removes duplicates based on read.id. The output file is written twice
    (once for concatenation, once after deduplication).
    """
    with open(fastq_out, 'w') as outfile:
        # concat
        with open(fastq_in1, 'r') as infile1:
            for record in SeqIO.parse(infile1, 'fastq'):
                SeqIO.write(record, outfile, 'fastq')

        with open(fastq_in2, 'r') as infile2:
            for record in SeqIO.parse(infile2, 'fastq'):
                SeqIO.write(record, outfile, 'fastq')

    dereplicate_fastq(fastq_in=fastq_out, fastq_out=fastq_out)


def get_chromosome(fasta: str, output_handle: str) -> None:
    """
    Extract the primary chromosome from a FASTA file.

    If the input contains a single contig, it is written to the output.
    If multiple contigs exist, the longest contig is selected and written
    as it is assumed to be the main chromosome. Logs information about
    the selected contig.

    Parameters
    ----------
    fasta : str
        Path to input FASTA file (single or multi-contig)
    output_handle : str
        Path for output FASTA file containing the selected chromosome

    Returns
    -------
    None
        Writes the selected chromosome to output_handle and logs the selection

    Notes
    -----
    Selection logic:
    - Single contig: Uses that contig directly
    - Multiple contigs: Selects the longest contig by sequence length

    The function assumes the longest contig is the main chromosome, which
    is appropriate for bacterial genomes or assemblies where the chromosome
    is expected to be significantly longer than plasmids or contaminants.

    Logging messages indicate which contig was selected and whether it was
    the only contig or chosen as the longest.
    """
    # test if there are a single entry in the fasta file
    try:  # there is a single entry
        chromosome = SeqIO.read(fasta, format='fasta')
        SeqIO.write(chromosome, output_handle, format='fasta')
        message = 'A single contig: {} was found and will be used for mapping'.format(
            '>' + chromosome.id
        )
        logging.info(message)

    except ValueError:  # there are more than one entry
        contigs = SeqIO.parse(fasta, format='fasta')
        max_len = 0

        # identify longest entry and assume it is the chromosome
        for record in contigs:
            seq_len = len(record.seq)

            if seq_len > max_len:
                chromosome = record
                max_len = seq_len

        SeqIO.write(chromosome, output_handle, format='fasta')
        message = 'The longest contig: {} has been saved as {} and will be used for mapping.'.format(
            '>' + chromosome.id, output_handle
        )
        logging.info(message)


def attach_seq(
    left: str, right: str, chromosome: str, output_name: str, offset: int = 0
) -> None:
    """
    Attach telomeric sequences to both ends of a chromosome sequence.

    Concatenates left and right sequences to the chromosome, optionally
    trimming bases from each end of the chromosome before attachment.
    This is used to build extended genomes with telomeric sequences.

    Parameters
    ----------
    left : str
        Path to FASTA file containing left/5' telomeric sequence
    right : str
        Path to FASTA file containing right/3' telomeric sequence
    chromosome : str
        Path to FASTA file containing chromosome sequence
    output_name : str
        Path for output FASTA file with attached sequences
    offset : int, default=0
        Number of bases to trim from each end of chromosome before attachment

    Returns
    -------
    None
        Writes extended genome to output_name

    Raises
    ------
    ValueError
        If offset is greater than or equal to half the chromosome length

    Notes
    -----
    The offset parameter allows trimming of chromosome ends to remove
    potentially problematic assembly regions before attaching telomeric
    sequences. If offset > 0, bases [offset:-offset] are retained.

    The output sequence ID is derived from output_name by removing the
    file extension.

    Example: For a 10kb chromosome with offset=100:
    - Chromosome bases 100-9900 are retained
    - Left sequence + chromosome[100:9900] + right sequence
    """
    left_seq = SeqIO.read(left, 'fasta')
    right_seq = SeqIO.read(right, 'fasta')
    chrom = SeqIO.read(chromosome, 'fasta')

    if offset == 0:  # if offset is 0 offset:-offset fucks it up
        genome = chrom
    elif offset >= len(chrom.seq) / 2:
        logging.error('Error: Offset is larger than  1/2 genome length.')
        return
    else:
        genome = chrom[offset:-offset]

    att_genome = left_seq + genome + right_seq
    att_genome.id = output_name.split('.')[0]
    SeqIO.write(att_genome, output_name, 'fasta')


# A function to merge fasta files
def merge_fasta(input_file1: str, input_file2: str, output_file: str) -> None:
    """
    Merge two FASTA files into a single multi-FASTA file.

    Combines all sequences from two FASTA files into one output file,
    preserving the order (file1 sequences first, then file2 sequences).
    Useful for creating multi-sequence reference files or combining
    consensus sequences.

    Parameters
    ----------
    input_file1 : str
        Path to first input FASTA file
    input_file2 : str
        Path to second input FASTA file
    output_file : str
        Path for output merged FASTA file

    Returns
    -------
    None
        Writes merged sequences to output_file

    Notes
    -----
    - All sequences from both files are included
    - Original sequence IDs and descriptions are preserved
    - Order is maintained: all sequences from file1, then all from file2
    - Can merge single-sequence or multi-sequence FASTA files
    """
    # Read sequences from input_file1 and input_file2
    sequences1 = list(SeqIO.parse(input_file1, 'fasta'))
    sequences2 = list(SeqIO.parse(input_file2, 'fasta'))

    # Merge sequences
    merged_sequences = sequences1 + sequences2

    # Write merged sequences to output_file
    with open(output_file, 'w') as output_handle:
        SeqIO.write(merged_sequences, output_handle, 'fasta')


def trim_to_cons(input_seq: str, num_base: int, output_handle: str) -> None:
    """
    Trim sequences to a specified number of bases from the start.

    Truncates all sequences in a FASTA file to the first num_base bases,
    adding a 'trimmed_' prefix to sequence IDs. Skips sequences shorter
    than the requested length with an error message.

    Parameters
    ----------
    input_seq : str
        Path to input FASTA file
    num_base : int
        Number of bases to retain from the start of each sequence
    output_handle : str
        Path for output trimmed FASTA file

    Returns
    -------
    None
        Writes trimmed sequences to output_handle

    Notes
    -----
    Processing details:
    - Sequences are trimmed to bases [0:num_base+1] (indices 0 through num_base)
    - Sequence IDs are prefixed with 'trimmed_'
    - Descriptions are removed from output sequences
    - Sequences shorter than num_base are skipped with error log
    - Only successfully trimmed sequences are written to output

    If all sequences are too short, an empty output file may be created.
    """
    # load file
    with open(input_seq) as fasta_file:
        all_rec = []

        for record in SeqIO.parse(fasta_file, 'fasta'):
            new_id = 'trimmed_' + record.id

            r_seq = record.seq
            length = len(r_seq)

            if num_base <= length:
                to_write = SeqRecord(
                    seq=r_seq[0 : num_base + 1], id=new_id, description=''
                )
                all_rec.append(to_write)

            else:
                logging.error('Error: Index out of range')

        if len(all_rec) > 0:
            SeqIO.write(all_rec, output_handle, 'fasta')


def strip_fasta(
    input_file: str, output_file: str, x: int, remove_from: str = 'start'
) -> None:
    """
    Remove a specified number of bases from sequence ends.

    Strips x bases from either the start (5' end) or end (3' end) of all
    sequences in a FASTA file. Useful for removing adapter sequences,
    low-quality ends, or trimming consensus sequences.

    Parameters
    ----------
    input_file : str
        Path to input FASTA file
    output_file : str
        Path for output stripped FASTA file
    x : int
        Number of bases to remove from each sequence
    remove_from : str, default='start'
        Which end to remove bases from: 'start' for 5' end, 'end' for 3' end

    Returns
    -------
    None
        Writes stripped sequences to output_file

    Raises
    ------
    AssertionError
        If x is not an integer
    ValueError
        If remove_from is not 'start' or 'end'

    Notes
    -----
    - Sequence IDs and descriptions are preserved
    - If remove_from='start': sequence[x:] is retained
    - If remove_from='end': sequence[:-x] is retained
    - All sequences in the file are processed identically
    - No validation that x is less than sequence length
    """
    assert type(x) is int

    records = []

    for record in SeqIO.parse(input_file, 'fasta'):
        if remove_from == 'start':
            modified_seq = record.seq[x:]
        elif remove_from == 'end':
            modified_seq = record.seq[:-x]
        else:
            raise ValueError("remove_from must be either 'start' or 'end'")

        record.seq = modified_seq
        records.append(record)

    SeqIO.write(records, output_file, 'fasta')


def build_extended_fasta(
    org_fasta: str, linear_elements: list[str], replicon_list: list, output_handle: str
) -> None:
    """
    Reconstruct multi-FASTA with extended linear contigs in original order.

    Replaces linear contigs that were extended by Telomore with their extended
    versions, while keeping circular/unprocessed contigs unchanged. The output
    maintains the original contig order and marks extended contigs as [linear].

    Parameters
    ----------
    org_fasta : str
        Path to original input FASTA file
    linear_elements : list of str
        List of contig IDs that were identified as linear and extended
    replicon_list : list
        List of Replicon objects containing paths to extended sequences
    output_handle : str
        Path for output FASTA file with extended contigs

    Returns
    -------
    None
        Writes reconstructed FASTA to output_handle

    Notes
    -----
    Processing logic:
    - Iterates through original FASTA in order
    - For linear contigs: replaces with extended version from Replicon.trim_out
    - For other contigs: copies unchanged from original
    - Adds '[linear]' to description of extended contigs

    This ensures the final assembly maintains the original contig order,
    which is important for tools that expect specific reference structures.
    The [linear] tag allows downstream tools to identify which contigs
    were extended.
    """
    seq_rec_list = []  # list of seqrecord to write to newfile

    for record in SeqIO.parse(org_fasta, 'fasta'):
        if record.id in linear_elements:
            for replicon in replicon_list:
                if replicon.name == record.id:
                    path_to_telomore_rec = replicon.trim_out
            telomore_rec = SeqIO.read(path_to_telomore_rec, format='fasta')
            telomore_rec.description = '[linear]'
            seq_rec_list.append(telomore_rec)
        else:
            seq_rec_list.append(record)

    SeqIO.write(sequences=seq_rec_list, handle=output_handle, format='fasta')


if __name__ == '__main__':
    check_fastq_order(
        '/tmp/tmpnrb8ke64/all_terminal_reads_1.fastq',
        '/tmp/tmpnrb8ke64/all_terminal_reads_2.fastq',
    )
