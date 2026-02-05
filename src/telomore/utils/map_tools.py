"""Functions for handling read mappings and extracting terminal reads."""

import gzip
import logging
from pathlib import Path
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam


def sam_to_readpair(
    sam_in: Path, fastq_in1: Path, fastq_in2: Path, fastq_out1: Path, fastq_out2: Path
) -> None:
    """
    Extract complete read pairs from paired-end FASTQ files based on SAM alignment.

    Retrieves both reads (R1 and R2) of paired-end sequences if either read
    appears in the input SAM file. This preserves read pairing for downstream
    analysis that requires synchronized paired-end data.

    Parameters
    ----------
    sam_in : Path
        Path to input SAM alignment file.
    fastq_in1 : Path
        Path to forward/R1 FASTQ file (gzip compressed).
    fastq_in2 : Path
        Path to reverse/R2 FASTQ file (gzip compressed).
    fastq_out1 : Path
        Path for output R1 FASTQ file containing extracted read pairs.
    fastq_out2 : Path
        Path for output R2 FASTQ file containing extracted read pairs.

    Returns
    -------
    None
        Writes extracted read pairs to fastq_out1 and fastq_out2.

    Notes
    -----
    Processing details:
    - First pass: Collects all read IDs from SAM file into a set
    - Handles read names with spaces by taking only the first part
    - Second pass: Extracts matching reads from gzipped FASTQ files
    - Both R1 and R2 are extracted if read name appears in SAM
    - Output files are uncompressed FASTQ format

    The input FASTQ files must be gzip-compressed (.gz), while outputs
    are plain text for immediate downstream processing.
    """
    with pysam.AlignmentFile(sam_in) as samfile:
        reads_to_grep = set()  # using a set should be faster than list

        # get all read names
        for read in samfile.fetch(until_eof=True):
            read_name = read.query_name
            if ' ' in read.query_name:
                read_name = read_name.split(' ')[0]

            reads_to_grep.add(read_name)

        # get read 1
        with (
            gzip.open(fastq_in1, 'rt') as gzip_handle,
            open(fastq_out1, 'w') as outfile,
        ):
            for record in SeqIO.parse(gzip_handle, 'fastq'):
                if record.id in reads_to_grep:
                    SeqIO.write(record, outfile, 'fastq')

        # get read 2
        with (
            gzip.open(fastq_in2, 'rt') as gzip_handle,
            open(fastq_out2, 'w') as outfile,
        ):
            for record in SeqIO.parse(gzip_handle, 'fastq'):
                if record.id in reads_to_grep:
                    SeqIO.write(record, outfile, 'fastq')


def sam_to_fastq(sam_in: Path, fastq_out: Path) -> None:
    r"""
    Convert SAM alignment file to FASTQ format, excluding unmapped reads.

    Extracts sequence and quality information from aligned reads in a SAM file
    and writes them in FASTQ format. Unmapped reads are filtered out. Used to
    extract reads that successfully aligned to terminal regions.

    Parameters
    ----------
    sam_in : Path
        Path to input SAM alignment file.
    fastq_out : Path
        File handle (opened in write mode) for output FASTQ.

    Returns
    -------
    None
        Writes FASTQ records to the provided file handle.

    Notes
    -----
    Processing details:
    - Only mapped reads (not flagged as unmapped) are converted
    - If quality scores are missing, assigns default high quality 'I' (Q40)
    - FASTQ format: @name\\nseq\\n+\\nqual\\n
    - The fastq_out parameter should be an open file handle, not a path string

    Quality score handling is important for reads extracted from SAM files
    that may not have retained the original quality information.
    """
    with pysam.AlignmentFile(sam_in, 'r') as samfile:
        for read in samfile.fetch(until_eof=True):
            if not read.is_unmapped:
                name = read.query_name
                seq = read.query_sequence
                qual = read.qual
                if qual is None:
                    qual = 'I' * len(seq)  # Assign a default high-quality score

                # Write the read in FASTQ format to the provided handle
                fastq_out.write(f'@{name}\n{seq}\n+\n{qual}\n')


def mapped_bases(cigarstring: str) -> int:
    """
    Calculate the number of bases mapped to the reference from a CIGAR string.

    Parses a CIGAR string and sums the lengths of all operations that consume
    reference bases (M, D, N, X, =). This is used to compare alignment quality
    when a read maps to multiple locations.

    Parameters
    ----------
    cigarstring : str
        CIGAR string from SAM alignment (e.g., '100M5S', '50M2D50M').

    Returns
    -------
    int
        Total number of bases that align to the reference sequence.

    Notes
    -----
    CIGAR operations that consume reference bases:
    - M: alignment match (can be match or mismatch)
    - D: deletion from reference
    - N: skipped region from reference
    - X: sequence mismatch
    - =: sequence match

    Operations that do NOT consume reference (excluded from count):
    - S: soft clipping
    - I: insertion to reference
    - H: hard clipping
    - P: padding

    This count represents how much of the reference sequence is covered
    by the alignment, which is useful for selecting the best alignment
    when a read maps to multiple positions.
    """
    # Define operations that consume reference bases
    consuming_operations = 'MDNX='

    # Parse the CIGAR string using regex
    # This produces a tuple in the format (121,"S")
    operations = re.findall(r'(\d+)([MIDNSHP=X])', cigarstring)

    # Initialize base count
    mapped_bases_count = 0

    # Loop through the parsed operations and sum bases for consuming operations
    for length, op in operations:
        if op in consuming_operations:
            mapped_bases_count += int(length)

    return mapped_bases_count


def cigar_maps_more_bases(cigar1: str, cigar2: str) -> bool:
    """
    Compare two CIGAR strings to determine which maps more reference bases.

    Evaluates which of two alignments covers more bases on the reference
    sequence. Used to select the better alignment when a read maps to
    multiple locations.

    Parameters
    ----------
    cigar1 : str
        First CIGAR string to compare.
    cigar2 : str
        Second CIGAR string to compare.

    Returns
    -------
    bool or None
        True if cigar1 maps more bases than cigar2, False if cigar2 maps
        more bases, None if they map equal bases.

    Notes
    -----
    The comparison is based on the number of reference-consuming bases
    (M, D, N, X, =) calculated by the mapped_bases function.

    Return values:
    - True: cigar1 has more mapped bases
    - False: cigar2 has more mapped bases
    - None: both have equal mapped bases (implicit, no return statement)

    This function is used to resolve multi-mapping reads by keeping the
    alignment that covers the most reference sequence, which typically
    indicates a better alignment quality.
    """
    bases1 = mapped_bases(cigar1)
    bases2 = mapped_bases(cigar2)

    if bases1 > bases2:
        return True
    elif bases1 < bases2:
        return False


def get_terminal_reads(
    sorted_bam_file: Path, contig: Path, loutput_handle: Path, routput_handle: Path
) -> None:
    """
    Extract reads mapping to the terminal 20bp regions of a contig.

    Retrieves all reads that align to the first or last 20 bases of a reference
    contig. For multi-mapping reads, keeps only the alignment with the most
    mapped bases. Critical for identifying reads that extend beyond assembly ends.

    Parameters
    ----------
    sorted_bam_file : Path
        Path to sorted BAM alignment file.
    contig : Path
        Name/ID of the contig to extract terminal reads from.
    loutput_handle : Path
        Path for output BAM file containing left-terminal reads.
    routput_handle : Path
        Path for output BAM file containing right-terminal reads.

    Returns
    -------
    None
        Writes left-terminal reads to loutput_handle and right-terminal reads
        to routput_handle.

    Notes
    -----
    Terminal region definition:
    - Left terminal: positions 0-20 (first 20bp)
    - Right terminal: (seq_end - 20) to seq_end (last 20bp)

    Multi-mapping read handling:
    - If a read maps to terminal region multiple times, compare CIGAR strings
    - Keep only the alignment mapping the most reference bases
    - Skips reads with no sequence (query_sequence is None)

    This function is essential for the Telomore workflow as it identifies
    reads that may contain sequence extending beyond the assembly, which
    can be used to build consensus extensions.
    """
    input = pysam.AlignmentFile(sorted_bam_file, 'r')

    # Fetch all reads aligned at start or end of reference
    seq_end = input.get_reference_length(contig)
    ref_name = contig
    left_reads = input.fetch(ref_name, start=0, stop=20)
    right_reads = input.fetch(ref_name, start=(seq_end - 20), stop=seq_end)

    # dict to store best mapped read from each end
    lterminal_reads = {}
    rterminal_reads = {}

    for lread in left_reads:
        query_name = lread.query_name
        cigar = lread.cigarstring

        if lread.query_sequence is None:  # skip empty reads
            continue

        # Check if the read is mapped multiple times and use
        # the read that maps to most bases
        if query_name in lterminal_reads:
            prior_read = lterminal_reads[query_name]
            prior_cigar = prior_read.cigarstring

            # Compare CIGAR strings to keep the one that maps more bases
            if cigar_maps_more_bases(cigar, prior_cigar):
                lterminal_reads[query_name] = lread
        else:
            lterminal_reads[query_name] = lread

    for rread in right_reads:
        query_name = rread.query_name
        cigar = rread.cigarstring

        if rread.query_sequence is None:  # skip empty reads
            continue

        # Check if the read is mapped multiple times and use
        # the read that maps to most bases
        if query_name in rterminal_reads:
            prior_read = rterminal_reads[query_name]
            prior_cigar = prior_read.cigarstring

            # Compare CIGAR strings to keep the one that maps more bases
            if cigar_maps_more_bases(cigar, prior_cigar):
                rterminal_reads[query_name] = rread
        else:
            rterminal_reads[query_name] = rread

    # Write all fetched reads to a new file
    lterminal_file = pysam.AlignmentFile(loutput_handle, 'w', template=input)
    for read in lterminal_reads.values():
        lterminal_file.write(read)
    lterminal_file.close()

    rterminal_file = pysam.AlignmentFile(routput_handle, 'w', template=input)
    for read in rterminal_reads.values():
        rterminal_file.write(read)
    rterminal_file.close()


def get_left_soft(sam_file: Path, left_out: Path, offset: int = 0) -> None:
    r"""
    Extract reads with 5' soft-clipping that extends beyond reference start.

    Identifies reads where the soft-clipped portion at the 5' end would extend
    beyond position 0 of the reference. Writes full alignments to SAM and
    extracts only the soft-clipped sequences to FASTQ. These represent sequence
    extending left of the assembly.

    Parameters
    ----------
    sam_file : Path
        Path to input SAM alignment file.
    left_out : Path
        Base path for output files (adds .sam and .fastq extensions).
    offset : int, default=0
        Additional bases to include beyond the soft-clipped region.

    Returns
    -------
    None
        Creates two output files:
        - {left_out}.sam: Full alignment records.
        - {left_out}.fastq: Soft-clipped sequences only.

    Notes
    -----
    Filtering logic:
    - Looks for CIGAR patterns starting with soft-clip: ^(\\d+)S
    - Only keeps reads where soft-clip length > reference_start position
    - This ensures the clipped sequence extends beyond the reference start

    FASTQ output contains:
    - Sequence: bases [0:clip_num+offset] from read
    - Quality: Phred scores converted to Sanger ASCII (Q+33)

    The offset parameter allows including additional bases for context,
    which can improve consensus building at the assembly boundary.
    """
    sam_in = pysam.AlignmentFile(sam_file, 'r')
    lclip = pysam.AlignmentFile(left_out + '.sam', 'w', template=sam_in)
    lfastq = open(left_out + '.fastq', 'w')

    start_clip = r'^(\d+)S'
    for read in sam_in:
        lmatch = re.match(start_clip, read.cigarstring)

        if lmatch:
            clip_num = int(lmatch.group(1))  # digits are retrieve via .group

            if clip_num > read.reference_start:
                lclip.write(read)  # write to sam-file

                # get info for fastq-file
                name = read.query_name
                seq = read.query_sequence[0 : clip_num + offset]
                sanger_qual = ''.join(
                    [chr(q + 33) for q in read.query_qualities[0 : clip_num + offset]]
                )  # phred qual converted to ASCII with 33 offset
                lfastq.write('@{}\n{}\n+\n{}\n'.format(name, seq, sanger_qual))
    sam_in.close()
    lclip.close()
    lfastq.close()


def get_right_soft(
    sam_file: Path, contig: Path, right_out: Path, offset: int = 0
) -> None:
    r"""
    Extract reads with 3' soft-clipping that extends beyond reference end.

    Identifies reads where the soft-clipped portion at the 3' end would extend
    beyond the reference sequence end. Writes full alignments to SAM and
    extracts only the soft-clipped sequences to FASTQ. These represent sequence
    extending right of the assembly.

    Parameters
    ----------
    sam_file : Path
        Path to input SAM alignment file.
    contig : Path
        Name/ID of the contig to determine reference length.
    right_out : Path
        Base path for output files (adds .sam and .fastq extensions).
    offset : int, default=0
        Additional bases to include beyond the soft-clipped region.

    Returns
    -------
    None
        Creates two output files:
        - {right_out}.sam: Full alignment records.
        - {right_out}.fastq: Soft-clipped sequences only.

    Notes
    -----
    Filtering logic:
    - Looks for CIGAR patterns ending with soft-clip: (\\d+)S$
    - Only keeps reads where (clip_length + reference_end) > seq_end
    - This ensures the clipped sequence extends beyond the reference end

    FASTQ output contains:
    - Sequence: last (clip_num+offset) bases from read
    - Quality: Phred scores converted to Sanger ASCII (Q+33)

    The offset parameter allows including additional bases for context,
    which can improve consensus building at the assembly boundary.
    """
    sam_in = pysam.AlignmentFile(sam_file, 'r')
    rclip = pysam.AlignmentFile(right_out + '.sam', 'w', template=sam_in)
    rfastq = open(right_out + '.fastq', 'w')
    seq_end = sam_in.get_reference_length(contig)  # get length of reference
    end_clip = r'(\d+)S$'
    for read in sam_in:
        rmatch = re.search(end_clip, read.cigarstring)
        if rmatch:
            clip_num = int(rmatch.group(1))  # digits are retrieve via .group

            if clip_num + read.reference_end > seq_end:
                rclip.write(read)  # write to sam-file

                # get info for fastq-file
                name = read.query_name
                seq = read.query_sequence[-(clip_num + offset) :]
                sanger_qual = ''.join(
                    [chr(q + 33) for q in read.query_qualities[-(clip_num + offset) :]]
                )  # phred qual converted to ASCII with 33 offset
                rfastq.write('@{}\n{}\n+\n{}\n'.format(name, seq, sanger_qual))

    sam_in.close()
    rclip.close()
    rfastq.close()


def revcomp_reads(reads_in: str, reads_out: str) -> None:
    """
    Generate reverse complement of all reads in a FASTQ file.

    Converts all sequences in a FASTQ file to their reverse complement,
    reversing both the sequence and quality scores. Adds 'rev_' prefix
    to read IDs. Used to orient left-terminal reads for consensus building.

    Parameters
    ----------
    reads_in : str
        Path to input FASTQ file.
    reads_out : str
        Path for output reverse-complemented FASTQ file.

    Returns
    -------
    None
        Writes reverse-complemented reads to reads_out.

    Notes
    -----
    Transformation details:
    - Sequence: Reverse complemented (A↔T, G↔C, reversed)
    - Quality scores: Reversed to match new sequence orientation
    - Read ID: Prefixed with 'rev_'
    - Original ID and quality annotations are preserved in structure

    This is necessary for left-terminal reads because they need to be
    reverse-complemented before consensus building to match the expected
    5' to 3' orientation for extension sequences.
    """
    with open(reads_in, 'r') as input_handle, open(reads_out, 'w') as output_handle:
        for record in SeqIO.parse(input_handle, 'fastq'):
            # Get the reverse complement of the sequence
            rev_complement_seq = record.seq.reverse_complement()

            # Reverse the quality scores as well
            rev_quality_scores = record.letter_annotations['phred_quality'][::-1]

            # Create a new record with the reverse complement sequence and quality scores
            rev_complement_record = record
            rev_complement_record.id = 'rev_' + str(record.id)

            rev_complement_record.seq = rev_complement_seq
            rev_complement_record.letter_annotations['phred_quality'] = (
                rev_quality_scores
            )

            # Write the reverse complement record to the output FASTQ file
            SeqIO.write(rev_complement_record, output_handle, 'fastq')


def revcomp(fasta_in: str, fasta_out: str) -> None:
    """
    Generate reverse complement of all sequences in a FASTA file.

    Converts all sequences in a FASTA file to their reverse complement.
    Adds 'rev_' prefix to sequence IDs. Used to reorient consensus sequences
    to match expected telomere orientation.

    Parameters
    ----------
    fasta_in : str
        Path to input FASTA file.
    fasta_out : str
        Path for output reverse-complemented FASTA file.

    Returns
    -------
    None
        Writes reverse-complemented sequences to fasta_out.

    Notes
    -----
    Transformation details:
    - Sequence: Reverse complemented (A↔T, G↔C, reversed)
    - Sequence ID: Prefixed with 'rev_'
    - Description preserved from original

    Unlike revcomp_reads, this operates on FASTA format and doesn't
    need to handle quality scores. Used primarily for consensus sequences
    built from left-terminal reads.
    """
    with open(fasta_in, 'r') as input_handle, open(fasta_out, 'w') as output_handle:
        for record in SeqIO.parse(input_handle, 'fasta'):
            # Get the reverse complement of the sequence
            rev_complement_seq = record.seq.reverse_complement()

            # Create a new record with the reverse complement sequence and quality scores
            rev_complement_record = record
            rev_complement_record.id = 'rev_' + str(record.id)
            rev_complement_record.seq = rev_complement_seq

            # Write the reverse complement record to the output FASTQ file
            SeqIO.write(rev_complement_record, output_handle, 'fasta')


def is_map_empty(file_path: str) -> bool:
    """
    Check if a BAM file contains any reads.

    Attempts to fetch the first read from a BAM alignment file to determine
    if the file is empty. Used to validate that alignment steps produced
    output before proceeding with downstream analysis.

    Parameters
    ----------
    file_path : str
        Path to BAM file to check.

    Returns
    -------
    bool
        False if the file contains at least one read, True if empty.

    Notes
    -----
    Implementation uses next() to attempt fetching the first read:
    - If successful: Returns False (file not empty)
    - If StopIteration raised: Returns True (file is empty)

    This is more efficient than loading all reads since it stops at
    the first read found. Empty BAM files indicate no reads aligned
    in a mapping step, which may require special handling.
    """
    # Open the alignment file
    with pysam.AlignmentFile(file_path, 'rb') as alignment_file:
        # Try to fetch the first read
        try:
            next(alignment_file)
            return False  # Alignment is not empty
        except StopIteration:
            return True  # Alignment is empty


def is_consensus_unmapped(file_path: str) -> bool:
    """
    Check if all reads in a BAM file are unmapped.

    Determines whether a consensus sequence failed to map to the reference
    by checking if all reads in the BAM file are flagged as unmapped. Used
    to detect when a consensus doesn't match the expected location.

    Parameters
    ----------
    file_path : str
        Path to BAM file to check.

    Returns
    -------
    bool
        True if all reads are unmapped or file is empty, False if any
        read is mapped.

    Notes
    -----
    Processing logic:
    - Loads all reads into memory (suitable for small consensus BAMs)
    - Returns True if file is empty (no reads)
    - Returns False immediately upon finding first mapped read
    - Returns True only if all reads are unmapped

    An unmapped consensus indicates the consensus sequence doesn't align
    to the expected position on the reference, suggesting it may not be
    a valid extension or may belong elsewhere in the genome.
    """
    with pysam.AlignmentFile(file_path, 'rb') as alignment_file:
        reads = list(alignment_file)  # get reads

        is_unmapped = True  # s

        if len(reads) > 0:
            for read in reads:
                if not read.is_unmapped:
                    is_unmapped = False
                    return is_unmapped

        return is_unmapped


def is_consensus_empty(file_path: str) -> bool:
    """
    Check if a BAM file represents an empty consensus sequence.

    Identifies BAM files produced by mapping empty consensus sequences, which
    contain exactly one unmapped read with no sequence. This indicates no
    consensus could be built, typically because no reads extended the assembly.

    Parameters
    ----------
    file_path : str
        Path to BAM file to check.

    Returns
    -------
    bool
        True if the file contains exactly one unmapped read with no sequence,
        False otherwise.

    Notes
    -----
    Criteria for empty consensus:
    1. Exactly one read in the file
    2. Read is flagged as unmapped
    3. Read has no sequence (seq is None or '*')

    This specific pattern occurs when an empty FASTA sequence (often produced
    when no terminal reads are found) is mapped against the reference. The
    aligner produces a single unmapped record with no sequence data.

    Distinguishes between:
    - Empty consensus: No reads to build consensus from
    - Unmapped consensus: Consensus built but doesn't align to expected location
    """
    with pysam.AlignmentFile(file_path, 'rb') as alignment_file:
        reads = list(alignment_file)  # Load all reads into a list

        # Check if there is exactly one read
        if len(reads) == 1:
            read = reads[0]
            # Check if the read is unmapped and has no sequence
            if read.is_unmapped and (not read.seq or read.seq == '*'):
                return True  # Only one unmapped read with no sequence
        return False  # Either more reads, or the read does not meet the conditions


def stitch_telo(
    ref: str,
    left_map: str,
    right_map: str,
    outfile: str,
    logout: str,
    tmp_left: str,
    tmp_right: str,
) -> tuple[int, int]:
    """
    Extend reference sequence with consensus sequences from terminal alignments.

    Extracts soft-clipped portions of consensus sequences that extend beyond
    the reference ends, attaches them to the reference, and creates a log
    documenting the extension process. Handles cases where consensus is empty,
    unmapped, or doesn't extend beyond reference.

    Parameters
    ----------
    ref : str
        Path to reference FASTA file.
    left_map : str
        Path to BAM file with left consensus aligned to reference.
    right_map : str
        Path to BAM file with right consensus aligned to reference.
    outfile : str
        Path for output extended FASTA file.
    logout : str
        Path for output log file documenting extension.
    tmp_left : str
        Path for temporary left consensus FASTA file.
    tmp_right : str
        Path for temporary right consensus FASTA file.

    Returns
    -------
    tuple of (int, int)
        Length of left consensus and length of right consensus in bases.

    Notes
    -----
    Left consensus processing:
    - Extracts reads mapping near reference start (position < 1000)
    - Looks for 5' soft-clipping extending beyond position 0
    - Adjusts for offset between soft-clip and actual overhang
    - Logs if consensus is empty, unmapped, or doesn't extend reference

    Right consensus processing:
    - Extracts reads mapping near reference end (position > 1000)
    - Looks for 3' soft-clipping extending beyond reference length
    - Adjusts for offset between soft-clip and actual overhang
    - Logs if consensus is empty, unmapped, or doesn't extend reference

    The output file contains: left_consensus + original_reference + right_consensus

    Empty SeqRecord objects are created when consensus fails validation,
    allowing the workflow to continue without breaking on concatenation.

    Log file format includes:
    - Section header
    - Consensus lengths
    - Error messages if consensus rejected
    - Full consensus sequences
    """
    left_log_mes = ''
    # Check if an empty left consensus was used to generate the map:
    if is_consensus_empty(left_map):
        # Make an empty seq list to enable errors later on
        left_seqs = []
        left_log_mes = '#No consensus produced for left-side end. Likely, no reads extends the assembly. '
    elif is_consensus_unmapped(left_map):
        left_seqs = []
        left_log_mes = f'#The consensus produced for the left-side does not map to left-side of {ref}'
    else:
        # extract left cons-to-stitch
        l_sam_in = pysam.AlignmentFile(left_map, 'r')
        left_seqs = []
        start_clip = r'^(\d+)S'
        # filter away mapping at right side
        cons_at_left = [read for read in l_sam_in if read.reference_start < 1000]

        # Get the sequence extending beyond the genome
        for read in cons_at_left:
            lmatch = re.match(start_clip, read.cigarstring)
            if lmatch:
                clip_num = int(lmatch.group(1))  # digits are retrieve via .group

                # check if the clipped sequence extends beyond genome
                if clip_num - read.reference_start <= 0:
                    left_log_mes = f'#The consensus produced for the left-side does extend beyond the start of {ref}'
                    left_seqs = []
                else:
                    seq = read.query_sequence[
                        0 : (clip_num - read.reference_start)
                    ]  # Adjust for if more than just overhanging bases are soft-clipped
                    left_seqs.append(seq)
        l_sam_in.close()

    right_log_mes = ''

    # Check if an empty left consensus was used to generate the map:
    if is_consensus_empty(right_map):
        right_seqs = []
        right_log_mes = '#No consensus produced for right-side end. Likely, no reads extends the assembly.'
    elif is_consensus_unmapped(right_map):
        right_seqs = []
        right_log_mes = f'#The consensus produced for the right-side does not map to the right-side of {ref}'
    else:
        # extract right cons-to-stitch
        r_sam_in = pysam.AlignmentFile(right_map, 'r')
        seq_end = r_sam_in.lengths[0]  # get length of reference
        right_seqs = []
        end_clip = r'(\d+)S$'  # reg. exp for ending with *S[num]

        cons_at_right = [read for read in r_sam_in if read.reference_start > 1000]
        for read in cons_at_right:
            rmatch = re.search(end_clip, read.cigarstring)
            if rmatch:
                clip_num = int(rmatch.group(1))  # digits are retrieve via .group
                # Adjusting for potential difference between overhang and soft-clip
                adj = seq_end - read.reference_end
                if clip_num + read.reference_end > seq_end:
                    seq = read.query_sequence[-(clip_num - adj) :]
                    right_seqs.append(seq)
        r_sam_in.close()

    # stitch the fuckers toghether
    genome = SeqIO.read(ref, 'fasta')

    # check if no conesnsus extens beyond the reference
    if len(left_seqs) == 0:
        left_cons = SeqRecord(
            Seq('')
        )  # if it is empty make an empty seqrecord to avoid errors in joining later
        logging.info('Left consensus does not extend genome')
    else:
        left_cons = SeqRecord(Seq(left_seqs[0]), id='left_cons')
        logging.info(f'Left consensus is {len(left_cons)}')
    if len(right_seqs) == 0:
        right_cons = SeqRecord(
            Seq('')
        )  # if it is empty make an empty seqrecord to avoid errors in joining later
        logging.info('Right cons does not extend genome')
    else:
        right_cons = SeqRecord(Seq(right_seqs[0]), id='right_cons')

        logging.info(f'Right consensus is {len(right_cons)}')
    new_genome = left_cons + genome + right_cons
    new_genome.id = 'Reference_with_consensus_attached'
    new_genome.description = ''
    SeqIO.write(new_genome, outfile, 'fasta')
    SeqIO.write(left_cons, tmp_left, 'fasta')
    SeqIO.write(right_cons, tmp_right, 'fasta')

    # Create log of consensus length
    log = open(logout, 'w')
    log.write(
        '=============================================================================='
    )
    log.write('\nINTIAL CONSENSUS')
    log.write(
        '\n=============================================================================='
    )
    log_content = '\nleft_cons:{}\tright_consensus:{}'.format(
        len(left_cons), len(right_cons)
    )
    comment_mes = '\n' + '\n'.join([left_log_mes, right_log_mes])
    log_content = log_content + comment_mes
    log.write(log_content)
    log.write('\n>left_cons\n')
    log.write(str(left_cons.seq))
    log.write('\n>right_cons\n')
    log.write(str(right_cons.seq))
    log.close()

    return (len(left_cons), len(right_cons))


def get_support_info(
    bam_file: str, genome: str, position: int, qual_threshold: int = 1
) -> tuple[int, int]:
    """
    Calculate coverage and reference-matching bases at a specific position.

    Determines read support at a genomic position by counting total coverage
    and the number of bases matching the reference. Used to validate consensus
    sequence quality by assessing read support at each position.

    Parameters
    ----------
    bam_file : str
        Path to BAM alignment file.
    genome : str
        Path to reference FASTA file.
    position : int
        Zero-based position to query.
    qual_threshold : int, default=1
        Minimum base quality score to include in counts.

    Returns
    -------
    tuple of (int, int)
        (coverage, matching_bases) where:
        - coverage: Total number of bases at this position.
        - matching_bases: Number of bases matching the reference.

    Notes
    -----
    Base counting:
    - Counts A, C, G, T bases separately at the position
    - Only includes bases with quality >= qual_threshold
    - Includes secondary mappings (read_callback='nofilter')
    - Sums all bases for total coverage

    Reference matching:
    - Compares reference base at position to read bases
    - If reference is 'N': matching_bases = 0
    - Otherwise: matching_bases = count of bases matching reference

    The matching ratio (matching_bases/coverage) indicates how well
    reads support the reference sequence at that position. High ratios
    (>0.7) indicate strong support, while low ratios suggest the consensus
    may not be well-supported by the reads.
    """
    fasta_file = SeqIO.read(genome, 'fasta')
    bam_in = pysam.AlignmentFile(
        bam_file,
        'rb',
    )

    # Set read_callback="no filter" to include secondary-mappings
    # Set quality threshold=1 to include all reads

    # Get reference name from BAM file
    reference_name = bam_in.get_reference_name(0)
    coverage_count = bam_in.count_coverage(
        reference_name,
        start=position,
        stop=position + 1,
        read_callback='nofilter',
        quality_threshold=qual_threshold,
    )
    A_num = coverage_count[0][0]
    C_num = coverage_count[1][0]
    G_num = coverage_count[2][0]
    T_num = coverage_count[3][0]
    cov = A_num + C_num + G_num + T_num

    if fasta_file.seq[position].upper() == 'N':
        matching_bases = 0
    elif fasta_file.seq[position].upper() == 'A':
        matching_bases = A_num
    elif fasta_file.seq[position].upper() == 'C':
        matching_bases = C_num
    elif fasta_file.seq[position].upper() == 'G':
        matching_bases = G_num
    elif fasta_file.seq[position].upper() == 'T':
        matching_bases = T_num

    return (cov, matching_bases)


def trim_by_map(
    untrimmed_assembly: str,
    sorted_bam_file: str,
    output_handle: str,
    cons_log: str,
    cov_thres: int = 5,
    ratio_thres: float = 0.7,
    qual_thres: int = 0,
) -> None:
    """
    Trim consensus extensions based on read support thresholds (Nanopore).

    Validates attached consensus sequences by trimming from the ends inward
    until finding positions with sufficient coverage and reference support.
    Removes unsupported consensus bases while retaining well-supported extensions.
    Optimized for Nanopore data with lower coverage requirements.

    Parameters
    ----------
    untrimmed_assembly : str
        Path to FASTA file with untrimmed consensus attached.
    sorted_bam_file : str
        Path to sorted BAM of terminal reads aligned to untrimmed assembly.
    output_handle : str
        Path for output trimmed FASTA file.
    cons_log : str
        Path to existing log file (will be appended with trimming info).
    cov_thres : int, default=5
        Minimum coverage depth required to keep a position.
    ratio_thres : float, default=0.7
        Minimum fraction of reads matching reference to keep a position.
    qual_thres : int, default=0
        Minimum base quality score to include in coverage calculation.

    Returns
    -------
    None
        Writes trimmed assembly to output_handle and appends to cons_log.

    Notes
    -----
    Trimming algorithm:
    1. Reads original consensus lengths from log file line 4
    2. Left end: Scans positions 0 to left_length
       - Stops at first position meeting coverage and ratio thresholds
       - Trims all bases before this position
    3. Right end: Scans positions (end - right_length) to end
       - Stops at first position meeting coverage and ratio thresholds
       - Trims all bases after this position

    Validation criteria:
    - Coverage >= cov_thres
    - (matching_bases / coverage) > ratio_thres
    - Base quality >= qual_thres

    Outcomes logged for each end:
    - Both rejected: Returns original reference only
    - One rejected: Keeps validated consensus on one side only
    - Both validated: Keeps both trimmed consensus sequences

    The output sequence ID indicates whether consensus was attached and
    includes descriptive suffix about trimming results.

    Designed for Nanopore data: Lower coverage threshold (5x) but
    similar ratio threshold to Illumina version.
    """
    # load genome
    fasta = SeqIO.read(untrimmed_assembly, 'fasta')
    fasta_end = len(fasta.seq) - 1  # subtract one to make it 0-indexed
    txt = open(cons_log, 'r')
    txt_lines = txt.readlines()[3]
    txt.close()
    left_len = int(txt_lines.split('\t')[0].split(':')[1])
    right_len = int(txt_lines.split('\t')[1].split(':')[1])

    index_start = None
    index_end = None

    # trim start/left-side
    for pos in range(0, 0 + left_len):
        try:
            cov, match = get_support_info(
                sorted_bam_file, untrimmed_assembly, pos, qual_thres
            )

            if cov >= cov_thres and (match / cov) > ratio_thres:
                index_start = pos

                break
        except TypeError:  # if no reads are mapped
            continue

    # trim end/right
    for pos in range(fasta_end, fasta_end - right_len, -1):
        try:
            cov, match = get_support_info(
                sorted_bam_file, untrimmed_assembly, pos, qual_thres
            )

            if cov >= cov_thres and (match / cov) > ratio_thres:
                index_end = pos

                break
        except TypeError:
            continue

    # check if coverage is too low for either consensus
    # Unclear on why, but adding one on the right side is nessesary to not trim an additional base
    # Even if the consensus is rejected.
    if index_start is None and index_end is None:
        trimmed_fasta = fasta[(0 + left_len) : (fasta_end - right_len) + 1]
        log_message = '\nLeft consensus rejected\nRight consensus rejected\n'
        trimmed_fasta.id = output_handle.split('.')[0] + '_with_no_consensus'
        trimmed_fasta.description = ''
    elif index_start is None:  # index without left consensus, but + right side
        log_message = (
            '\nLeft consensus rejected\nRight consensus trimmed with {}\n'.format(
                (fasta_end - index_end)
            )
        )
        trimmed_fasta = fasta[(0 + left_len) : index_end + 1]
        trimmed_fasta.id = (
            output_handle.split('.')[0] + '_with_trimmed_consensus_attached'
        )
        trimmed_fasta.description = ''
    elif index_end is None:  # index from consensus until before consensus on right side
        log_message = '\nLeft consensus trimmed with {}\nRight rejected\n'.format(
            index_start
        )
        trimmed_fasta = fasta[index_start : (fasta_end - right_len) + 1]
        trimmed_fasta.id = (
            output_handle.split('.')[0] + '_with_trimmed_consensus_attached'
        )
        trimmed_fasta.description = ''
    else:
        log_message = '\nLeft consensus trimmed with {}\nRight consensus trimmed with {}\n'.format(
            index_start, (fasta_end - index_end)
        )
        trimmed_fasta = fasta[index_start : index_end + 1]
        trimmed_fasta.id = (
            output_handle.split('.')[0] + '_with_trimmed_consensus_attached'
        )
        trimmed_fasta.description = ''

    log = open(cons_log, 'a')
    log.write(
        '\n=============================================================================='
    )
    log.write('\nCONSENSUS TRIMMING')
    log.write(
        '\n=============================================================================='
    )
    log.write(
        f'\nRule: Trimmed until Q_score>= {qual_thres}, cov>= {cov_thres} and supporting ratio>= {ratio_thres}'
    )
    log.write(log_message)
    log.close()
    SeqIO.write(trimmed_fasta, output_handle, 'fasta')


def trim_by_map_illumina(
    untrimmed_assembly: str,
    sorted_bam_file: str,
    output_handle: str,
    cons_log: str,
    cov_thres: int = 1,
    ratio_thres: float = 0.7,
    qual_thres: int = 30,
) -> None:
    """
    Trim consensus extensions based on read support thresholds (Illumina).

    Validates attached consensus sequences by trimming from the ends inward
    until finding positions with sufficient coverage and reference support.
    Removes unsupported consensus bases while retaining well-supported extensions.
    Optimized for Illumina data with high quality requirements.

    Parameters
    ----------
    untrimmed_assembly : str
        Path to FASTA file with untrimmed consensus attached.
    sorted_bam_file : str
        Path to sorted BAM of terminal reads aligned to untrimmed assembly.
    output_handle : str
        Path for output trimmed FASTA file.
    cons_log : str
        Path to existing log file (will be appended with trimming info).
    cov_thres : int, default=1
        Minimum coverage depth required to keep a position.
    ratio_thres : float, default=0.7
        Minimum fraction of reads matching reference to keep a position.
    qual_thres : int, default=30
        Minimum base quality score (Q30) to include in coverage calculation.

    Returns
    -------
    None
        Writes trimmed assembly to output_handle and appends to cons_log.

    Notes
    -----
    Trimming algorithm:
    1. Reads original consensus lengths from log file line 4
    2. Left end: Scans positions 0 to left_length
       - Stops at first position meeting coverage and ratio thresholds
       - Trims all bases before this position
    3. Right end: Scans positions (end - right_length) to end
       - Stops at first position meeting coverage and ratio thresholds
       - Trims all bases after this position

    Validation criteria:
    - Coverage >= cov_thres
    - (matching_bases / coverage) > ratio_thres
    - Base quality >= qual_thres

    Outcomes logged for each end:
    - Both rejected: Returns original reference only
    - One rejected: Keeps validated consensus on one side only
    - Both validated: Keeps both trimmed consensus sequences

    The output sequence ID indicates whether consensus was attached and
    includes descriptive suffix about trimming results.

    Designed for Illumina data: Higher quality threshold (Q30) but
    lower coverage requirement (1x) compared to Nanopore version.
    Illumina's higher per-base accuracy allows more stringent quality
    filtering with lower coverage depth.
    """
    # load genome
    fasta = SeqIO.read(untrimmed_assembly, 'fasta')
    fasta_end = len(fasta.seq) - 1  # subtract one to make it 0-indexed
    txt = open(cons_log, 'r')
    txt_lines = txt.readlines()[3]
    txt.close()
    left_len = int(txt_lines.split('\t')[0].split(':')[1])
    right_len = int(txt_lines.split('\t')[1].split(':')[1])

    index_start = None
    index_end = None

    # trim start/left-side
    for pos in range(0, 0 + left_len):
        try:
            cov, match = get_support_info(
                sorted_bam_file, untrimmed_assembly, pos, qual_thres
            )

            if cov >= cov_thres and (match / cov) > ratio_thres:
                index_start = pos

                break
        except TypeError:  # if no reads are mapped
            continue

    # trim end/right
    for pos in range(fasta_end, fasta_end - right_len, -1):
        try:
            cov, match = get_support_info(
                sorted_bam_file, untrimmed_assembly, pos, qual_thres
            )

            if cov >= cov_thres and (match / cov) > ratio_thres:
                index_end = pos

                break
        except TypeError:
            continue

    # check if coverage is too low for either consensus
    # Unclear on why, but adding one on the right side is nessesary to not trim an additional base
    # Even if the consensus is rejected.
    if index_start is None and index_end is None:
        trimmed_fasta = fasta[(0 + left_len) : (fasta_end - right_len) + 1]
        log_message = '\nLeft consensus rejected\nRight consensus rejected\n'
        trimmed_fasta.id = output_handle.split('.')[0] + '_with_no_consensus'
        trimmed_fasta.description = ''
    elif index_start is None:  # index without left consensus, but + right side
        log_message = (
            '\nLeft consensus rejected\nRight consensus trimmed with {}\n'.format(
                (fasta_end - index_end)
            )
        )
        trimmed_fasta = fasta[(0 + left_len) : index_end + 1]
        trimmed_fasta.id = (
            output_handle.split('.')[0] + '_with_trimmed_consensus_attached'
        )
        trimmed_fasta.description = ''
    elif index_end is None:  # index from consensus until before consensus on right side
        log_message = '\nLeft consensus trimmed with {}\nRight rejected\n'.format(
            index_start
        )
        trimmed_fasta = fasta[index_start : (fasta_end - right_len) + 1]
        trimmed_fasta.id = (
            output_handle.split('.')[0] + '_with_trimmed_consensus_attached'
        )
        trimmed_fasta.description = ''
    else:
        log_message = '\nLeft consensus trimmed with {}\nRight consensus trimmed with {}\n'.format(
            index_start, (fasta_end - index_end)
        )
        trimmed_fasta = fasta[index_start : index_end + 1]
        trimmed_fasta.id = (
            output_handle.split('.')[0] + '_with_trimmed_consensus_attached'
        )
        trimmed_fasta.description = ''

    log = open(cons_log, 'a')
    log.write(
        '\n=============================================================================='
    )
    log.write('\nCONSENSUS TRIMMING')
    log.write(
        '\n=============================================================================='
    )
    log.write(
        f'\nRule: Trimmed until Q_score>= {qual_thres}, cov>= {cov_thres} and supporting ratio>= {ratio_thres}'
    )
    log.write(log_message)
    log.close()
    SeqIO.write(trimmed_fasta, output_handle, 'fasta')


def generate_support_log(genome: str, qc_bam_file: str, output_handle: str) -> None:
    """
    Generate position-by-position coverage and support statistics for QC.

    Creates a detailed log showing coverage and reference-matching bases at
    every position in the genome. Used for quality control visualization and
    analysis of read support across the extended assembly.

    Parameters
    ----------
    genome : str
        Path to reference genome FASTA file.
    qc_bam_file : str
        Path to BAM file with QC reads aligned to genome.
    output_handle : str
        Path for output log file with coverage statistics.

    Returns
    -------
    None
        Writes position, coverage, and matching bases to output_handle.

    Notes
    -----
    For each position from 0 to (genome_length - 1):
    - Calculates coverage (total bases)
    - Calculates matching bases (bases matching reference)
    - Prints position, coverage, matching_bases to stdout
    - Skips positions where no reads map (TypeError caught)

    The output allows plotting coverage profiles to visualize:
    - Read support across the genome
    - Quality of consensus extensions at telomeres
    - Positions where support drops (potential trimming sites)

    Uses qual_threshold=1 to include all bases regardless of quality,
    providing a complete picture of coverage for QC purposes.

    Note: Current implementation only prints to stdout. To write to file,
    the log.write() call should be corrected.
    """
    # trim start/left-side

    fasta = SeqIO.read(genome, 'fasta')
    fasta_end = len(fasta.seq) - 1  # subtract one to make it 0-indexed

    # Generate log of coverage at all positions
    with open(output_handle, 'a') as log:
        for pos in range(0, fasta_end):
            try:
                cov, match = get_support_info(
                    bam_file=qc_bam_file, genome=genome, position=pos, qual_threshold=1
                )

                print(pos, cov, match)
                log.write(pos, cov, match)
            except TypeError:  # if no reads are mapped
                continue
