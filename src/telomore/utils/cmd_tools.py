"""
Functions for running CLI-tools to map reads and generate consensus.
"""

import logging
import os
from pathlib import Path
import subprocess
import traceback

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def map_and_sort(reference: str, fastq: str, output: str, threads: int = 1) -> None:
    """
    Map Oxford Nanopore long reads to a reference genome using minimap2.

    This function performs the complete mapping workflow for long reads:
    maps reads using minimap2 with parameters optimized for telomere extension
    analysis, sorts the output to BAM format, and creates an index file.

    Parameters
    ----------
    reference : str
        Path to the reference genome file (.fasta, .fna, or .fa)
    fastq : str
        Path to the Oxford Nanopore FASTQ file (can be gzipped)
    output : str
        Path for the output BAM file
    threads : int, default=1
        Number of threads to use for mapping and sorting

    Returns
    -------
    None
        Creates a sorted, indexed BAM file at the specified output path

    Raises
    ------
    FileNotFoundError
        If reference or FASTQ files don't exist
    subprocess.CalledProcessError
        If any of the mapping steps fail
    TypeError
        If threads is not an integer

    Notes
    -----
    Uses minimap2 with these key parameters for long read mapping:
    - `-a`: Output in SAM format
    - `--secondary-seq`: Include sequence in secondary alignments
    - `-Y`: Use soft clipping for supplementary alignments

    These parameters are crucial for telomere extension analysis as they ensure:
    - Secondary alignments retain sequence information
    - Soft clips are preserved for reads extending beyond genome edges

    The function automatically sorts the output and creates a BAM index.
    """
    # Input validation
    if not isinstance(threads, int):
        raise TypeError('threads must be an integer')
    if not os.path.isfile(reference):
        raise FileNotFoundError(f'Reference file does not exist: {reference}')
    if not os.path.isfile(fastq):
        raise FileNotFoundError(f'FASTQ file does not exist: {fastq}')

    # Define path for temporary SAM file
    temp_sam = f'{output}.temp.sam'

    logging.info(f'Starting Oxford Nanopore read mapping for {fastq}')
    logging.info(f'Reference: {reference}, Output: {output}, Threads: {threads}')

    try:
        # Step 1: Map reads using minimap2
        logging.info('Mapping long reads with minimap2...')
        map_cmd = [
            'minimap2',
            '-a',  # Output in SAM format
            reference,  # Reference genome
            fastq,  # Input FASTQ file
            '-t',
            str(threads),  # Number of threads
            '-o',
            temp_sam,  # Output SAM file
            '--secondary-seq',  # Include sequence in secondary alignments
            '-Y',  # Use soft clipping for supplementary alignments
        ]

        result = subprocess.run(map_cmd, capture_output=True, text=True, check=True)
        logging.debug(f'Minimap2 mapping completed: {result.stderr}')

        # Step 2: Sort SAM to BAM
        logging.info('Sorting SAM file to BAM format...')
        sort_cmd = [
            'samtools',
            'sort',
            '-@',
            str(threads),  # Number of threads
            '-o',
            output,  # Output file
            temp_sam,  # Input SAM file
        ]

        result = subprocess.run(sort_cmd, capture_output=True, text=True, check=True)
        logging.debug(f'SAM sorting completed: {result.stderr}')

        # Step 3: Index the BAM file
        logging.info('Indexing BAM file...')
        index_cmd = ['samtools', 'index', output]

        result = subprocess.run(index_cmd, capture_output=True, text=True, check=True)
        logging.debug(f'BAM indexing completed: {result.stderr}')

        logging.info(f'Oxford Nanopore mapping completed successfully: {output}')

    except subprocess.CalledProcessError as e:
        logging.error(f'Nanopore mapping failed at command: {" ".join(e.cmd)}')
        logging.error(f'Error code: {e.returncode}')
        logging.error(f'Stdout: {e.stdout}')
        logging.error(f'Stderr: {e.stderr}')
        raise

    except Exception as e:
        logging.error(f'Unexpected error during Nanopore mapping: {e}')
        logging.error(traceback.format_exc())
        raise

    finally:
        # Clean up temporary SAM file
        logging.info('Cleaning up temporary files...')
        if os.path.exists(temp_sam):
            try:
                os.remove(temp_sam)
                logging.debug(f'Removed temporary SAM file: {temp_sam}')
            except OSError as e:
                logging.warning(f'Could not remove temporary SAM file {temp_sam}: {e}')


def map_and_sort_illumina(
    reference: str, read1: str, read2: str, output: str, threads: int = 1
) -> None:
    """
    Map Illumina paired-end reads to a reference genome using Bowtie2.

    This function performs the complete mapping workflow: builds a Bowtie2 index,
    maps reads with parameters optimized for telomere extension analysis,
    sorts the output to BAM format, and creates an index file.

    Parameters
    ----------
    reference : str
        Path to the reference genome file (.fasta, .fna, or .fa)
    read1 : str
        Path to the first paired-end FASTQ file (gzipped)
    read2 : str
        Path to the second paired-end FASTQ file (gzipped)
    output : str
        Path for the output BAM file
    threads : int, default=1
        Number of threads to use for mapping and sorting

    Returns
    -------
    None
        Creates a sorted, indexed BAM file at the specified output path

    Raises
    ------
    FileNotFoundError
        If reference or read files don't exist
    subprocess.CalledProcessError
        If any of the mapping steps fail
    TypeError
        If threads is not an integer

    Notes
    -----
    Uses Bowtie2 with these key parameters:
    - `-a`: Report all alignments (including secondary/supplementary)
    - `--local`: Allow soft-clipping for reads extending beyond genome edges
    - `--no-mixed`: Suppress unpaired alignments for paired reads
    - `--no-discordant`: Suppress discordant alignments

    Temporary index files are automatically cleaned up after mapping.
    """
    # Input validation
    if not isinstance(threads, int):
        raise TypeError('threads must be an integer')
    if not os.path.isfile(reference):
        raise FileNotFoundError(f'Reference file does not exist: {reference}')
    if not os.path.isfile(read1):
        raise FileNotFoundError(f'Read1 file does not exist: {read1}')
    if not os.path.isfile(read2):
        raise FileNotFoundError(f'Read2 file does not exist: {read2}')

    # Define paths for temporary files
    ref_path = Path(reference)
    index_prefix = f'{ref_path}.bt.index'
    temp_sam = f'{output}.temp.sam'

    logging.info(f'Starting Illumina read mapping for {read1} and {read2}')
    logging.info(f'Reference: {reference}, Output: {output}, Threads: {threads}')

    try:
        # Step 1: Build Bowtie2 index
        logging.info('Building Bowtie2 index...')
        build_cmd = [
            'bowtie2-build',
            '-q',  # Quiet mode
            '--threads',
            str(threads),
            reference,
            index_prefix,
        ]

        result = subprocess.run(build_cmd, capture_output=True, text=True, check=True)
        logging.debug(f'Index build completed: {result.stderr}')

        # Step 2: Map reads using Bowtie2
        logging.info('Mapping reads with Bowtie2...')
        map_cmd = [
            'bowtie2',
            '-a',  # Report all alignments
            '-p',
            str(threads),  # Number of threads
            '--local',  # Local alignment mode for soft-clipping
            '--sam-no-qname-trunc',  # Don't truncate QNAME
            '-x',
            index_prefix,  # Index prefix
            '-1',
            read1,  # First mate file
            '-2',
            read2,  # Second mate file
            '-S',
            temp_sam,  # Output SAM file
            '--quiet',  # Suppress verbose output
            '--no-mixed',  # Suppress unpaired alignments
            '--no-discordant',  # Suppress discordant alignments
        ]

        result = subprocess.run(map_cmd, capture_output=True, text=True, check=True)
        logging.debug(f'Read mapping completed: {result.stderr}')

        # Step 3: Sort SAM to BAM
        logging.info('Sorting SAM file to BAM format...')
        sort_cmd = [
            'samtools',
            'sort',
            '-@',
            str(threads),  # Number of threads
            '-o',
            output,  # Output file
            temp_sam,  # Input SAM file
        ]

        result = subprocess.run(sort_cmd, capture_output=True, text=True, check=True)
        logging.debug(f'SAM sorting completed: {result.stderr}')

        # Step 4: Index the BAM file
        logging.info('Indexing BAM file...')
        index_cmd = ['samtools', 'index', output]

        result = subprocess.run(index_cmd, capture_output=True, text=True, check=True)
        logging.debug(f'BAM indexing completed: {result.stderr}')

        logging.info(f'Illumina mapping completed successfully: {output}')

    except subprocess.CalledProcessError as e:
        logging.error(f'Illumina mapping failed at command: {" ".join(e.cmd)}')
        logging.error(f'Error code: {e.returncode}')
        logging.error(f'Stdout: {e.stdout}')
        logging.error(f'Stderr: {e.stderr}')
        raise

    except Exception as e:
        logging.error(f'Unexpected error during Illumina mapping: {e}')
        logging.error(traceback.format_exc())
        raise

    finally:
        # Clean up temporary files
        logging.info('Cleaning up temporary files...')

        # Remove temporary SAM file
        if os.path.exists(temp_sam):
            try:
                os.remove(temp_sam)
                logging.debug(f'Removed temporary SAM file: {temp_sam}')
            except OSError as e:
                logging.warning(f'Could not remove temporary SAM file {temp_sam}: {e}')

        # Remove Bowtie2 index files
        index_extensions = [
            '.1.bt2',
            '.2.bt2',
            '.3.bt2',
            '.4.bt2',
            '.rev.1.bt2',
            '.rev.2.bt2',
        ]
        for ext in index_extensions:
            index_file = f'{index_prefix}{ext}'
            if os.path.exists(index_file):
                try:
                    os.remove(index_file)
                    logging.debug(f'Removed index file: {index_file}')
                except OSError as e:
                    logging.warning(f'Could not remove index file {index_file}: {e}')


def map_and_sort_illumina_cons(
    reference: str, consensus_fasta: str, output: str, threads: int = 1
) -> None:
    """
    Map consensus sequences to a reference genome using Bowtie2.

    This function performs the complete mapping workflow for consensus sequences:
    builds a Bowtie2 index, maps the consensus FASTA sequences with parameters
    optimized for telomere extension analysis, sorts the output to BAM format,
    and creates an index file.

    Parameters
    ----------
    reference : str
        Path to the reference genome file (.fasta, .fna, or .fa)
    consensus_fasta : str
        Path to the consensus sequences file in FASTA format
    output : str
        Path for the output BAM file
    threads : int, default=1
        Number of threads to use for mapping and sorting

    Returns
    -------
    None
        Creates a sorted, indexed BAM file at the specified output path

    Raises
    ------
    FileNotFoundError
        If reference or consensus files don't exist
    subprocess.CalledProcessError
        If any of the mapping steps fail
    TypeError
        If threads is not an integer

    Notes
    -----
    Uses Bowtie2 with these key parameters for consensus mapping:
    - `-a`: Report all alignments (including secondary/supplementary)
    - `--local`: Allow soft-clipping for sequences extending beyond genome edges
    - `-f`: Input sequences are in FASTA format (not FASTQ)
    - `--sam-no-qname-trunc`: Don't truncate sequence names

    This function is specifically designed for mapping consensus sequences
    generated from read assemblies, not raw sequencing reads.

    Temporary index files are automatically cleaned up after mapping.
    """
    # Input validation
    if not isinstance(threads, int):
        raise TypeError('threads must be an integer')
    if not os.path.isfile(reference):
        raise FileNotFoundError(f'Reference file does not exist: {reference}')
    if not os.path.isfile(consensus_fasta):
        raise FileNotFoundError(
            f'Consensus FASTA file does not exist: {consensus_fasta}'
        )

    # Define paths for temporary files
    ref_path = Path(reference)
    index_prefix = f'{ref_path}.bt.cons.index'
    temp_sam = f'{output}.temp.sam'

    logging.info(f'Starting consensus sequence mapping for {consensus_fasta}')
    logging.info(f'Reference: {reference}, Output: {output}, Threads: {threads}')

    try:
        # Step 1: Build Bowtie2 index
        logging.info('Building Bowtie2 index for consensus mapping...')
        build_cmd = [
            'bowtie2-build',
            '-q',  # Quiet mode
            '--threads',
            str(threads),
            reference,
            index_prefix,
        ]

        result = subprocess.run(build_cmd, capture_output=True, text=True, check=True)
        logging.debug(f'Index build completed: {result.stderr}')

        # Step 2: Map consensus sequences using Bowtie2
        logging.info('Mapping consensus sequences with Bowtie2...')
        map_cmd = [
            'bowtie2',
            '-a',  # Report all alignments
            '-p',
            str(threads),  # Number of threads
            '--local',  # Local alignment mode for soft-clipping
            '--sam-no-qname-trunc',  # Don't truncate QNAME
            '-f',  # Input is FASTA format (not FASTQ)
            '-x',
            index_prefix,  # Index prefix
            '-U',
            consensus_fasta,  # Unpaired input file (consensus sequences)
            '-S',
            temp_sam,  # Output SAM file
            '--quiet',  # Suppress verbose output
        ]

        result = subprocess.run(map_cmd, capture_output=True, text=True, check=True)
        logging.debug(f'Consensus mapping completed: {result.stderr}')

        # Step 3: Sort SAM to BAM
        logging.info('Sorting SAM file to BAM format...')
        sort_cmd = [
            'samtools',
            'sort',
            '-@',
            str(threads),  # Number of threads
            '-o',
            output,  # Output file
            temp_sam,  # Input SAM file
        ]

        result = subprocess.run(sort_cmd, capture_output=True, text=True, check=True)
        logging.debug(f'SAM sorting completed: {result.stderr}')

        # Step 4: Index the BAM file
        logging.info('Indexing BAM file...')
        index_cmd = ['samtools', 'index', output]

        result = subprocess.run(index_cmd, capture_output=True, text=True, check=True)
        logging.debug(f'BAM indexing completed: {result.stderr}')

        logging.info(f'Consensus mapping completed successfully: {output}')

    except subprocess.CalledProcessError as e:
        logging.error(f'Consensus mapping failed at command: {" ".join(e.cmd)}')
        logging.error(f'Error code: {e.returncode}')
        logging.error(f'Stdout: {e.stdout}')
        logging.error(f'Stderr: {e.stderr}')
        raise

    except Exception as e:
        logging.error(f'Unexpected error during consensus mapping: {e}')
        logging.error(traceback.format_exc())
        raise

    finally:
        # Clean up temporary files
        logging.info('Cleaning up temporary files...')

        # Remove temporary SAM file
        if os.path.exists(temp_sam):
            try:
                os.remove(temp_sam)
                logging.debug(f'Removed temporary SAM file: {temp_sam}')
            except OSError as e:
                logging.warning(f'Could not remove temporary SAM file {temp_sam}: {e}')

        # Remove Bowtie2 index files
        index_extensions = [
            '.1.bt2',
            '.2.bt2',
            '.3.bt2',
            '.4.bt2',
            '.rev.1.bt2',
            '.rev.2.bt2',
        ]
        for ext in index_extensions:
            index_file = f'{index_prefix}{ext}'
            if os.path.exists(index_file):
                try:
                    os.remove(index_file)
                    logging.debug(f'Removed index file: {index_file}')
                except OSError as e:
                    logging.warning(f'Could not remove index file {index_file}: {e}')


def train_lastDB(fasta_name: str, reads: str, db_name: str, t=1) -> None:
    """Trains and lastDB database using a reference and long-reads"""
    # index fasta file
    try:
        subprocess.run(
            ['lastdb', '-P' + str(t), '-uRY4', db_name, fasta_name], check=True
        )
    except subprocess.CalledProcessError as e:
        # If the bash script fails, capture the error and log the traceback
        logging.error('train_lastDB failed with error: %s', e)
        logging.error('Script stderr: %s', e.stderr)
        logging.error(traceback.format_exc())
        raise
    # train last-db
    try:
        file = open(db_name + '.par', 'w')  # needed to write to file
        subprocess.run(
            ['last-train', '-P' + str(t), '-Qkeep', db_name, reads],
            stdout=file,
            check=True,
        )
        file.close()
    except subprocess.CalledProcessError as e:
        # If the bash script fails, capture the error and log the traceback
        logging.error('train_lastDB failed with error: %s', e)
        logging.error('Script stderr: %s', e.stderr)
        logging.error(traceback.format_exc())
        raise


def generate_consensus_lamassemble(db_name: str, reads: str, output: str) -> None:
    """Generates a consensus fasta-file given a LAST-DB and a series of reads.
    Is suitable for dissimilar reads such as nanopore"""

    # Check if only a single read is present before generating consensus
    # If that is the case write just that read to consensus file
    sequence_count = 0
    for record in SeqIO.parse(reads, 'fastq'):
        sequence_count += 1
        latest_record = record

    if sequence_count == 0:
        logging.info(
            'There are no reads to construct a consensus from. Emtpy consensus returned to %s',
            output,
        )
        with open(f'{output}', 'w') as seq:
            empty_record = SeqRecord(Seq(''), id='empty_consensus')
            SeqIO.write(empty_record, seq, 'fasta')
    if sequence_count == 1:
        single_record = latest_record
        logging.info(
            'There are only a single read to construct a consensus from. Returning read as consensus to %s',
            output,
        )
        with open(f'{output}', 'w') as seq:
            SeqIO.write(single_record, seq, 'fasta')
    elif sequence_count > 1:
        db = db_name + '.par'
        seq = open(str(output), 'w')
        try:
            subprocess.run(
                ['lamassemble', '--name=' + output, db, reads], stdout=seq, check=True
            )
            seq.close()

            aln = open(str(output) + '.aln', 'w')  # save alignment of reads
            subprocess.run(['lamassemble', '-a', db, reads], stdout=aln)
            aln.close()

        except subprocess.CalledProcessError as e:
            # If the bash script fails, capture the error and log the traceback
            logging.error('generate_consensus_lamassemble failed with error: %s', e)
            logging.error('Script stderr: %s', e.stderr)
            logging.error(traceback.format_exc())
            raise


def generate_consensus_mafft(reads: str, output: str) -> None:
    """Generates a consensus fasta-file given a series of reads in fasta-format.
    Relies on very similar reads."""

    # Check if only a single read is mapped and use that as the consensus if there are no others.
    sequence_count = 0
    for record in SeqIO.parse(reads, 'fastq'):
        sequence_count += 1
        latest_record = record

    if sequence_count == 0:
        logging.info(
            'There are no reads to construct a consensus from. Emtpy consensus returned to %s',
            output,
        )
        with open(f'{output}', 'w') as seq:
            empty_record = SeqRecord(Seq(''), id='empty_consensus')
            SeqIO.write(empty_record, seq, 'fasta')

    if sequence_count == 1:
        single_record = latest_record
        logging.info(
            'There are only a single read to construct a consensus from. Returning read as consensus to %s',
            output,
        )
        with open(f'{output}', 'w') as seq:
            SeqIO.write(single_record, seq, 'fasta')

    elif sequence_count > 1:
        # Convert fastq to fasta
        fasta_reads = reads + '.fasta'
        with open(reads, 'r') as input_handle, open(fasta_reads, 'w') as output_handle:
            SeqIO.convert(input_handle, 'fastq', output_handle, 'fasta')
        # Run Mafft
        mafft_output = output + '.aln'
        try:
            mafft_file = open(mafft_output, 'w')
            subprocess.run(
                ['mafft', '--quiet', fasta_reads], stdout=mafft_file, check=True
            )
            mafft_file.close()
            os.remove(fasta_reads)

        except subprocess.CalledProcessError as e:
            # If the bash script fails, capture the error and log the traceback
            logging.error('generate_consensus_mafft failed with error: %s', e)
            logging.error('Script stderr: %s', e.stderr)
            logging.error(traceback.format_exc())
            raise

        try:
            # Generate consensus using Emboss cons
            # Plurality 1 ensures that only one reads needs to cover a position to generate
            # consensus. This is dangerous with dissimilar seqeunces
            subprocess.run(
                [
                    'cons',
                    '-name=' + output,
                    '-plurality',
                    str(1),
                    '-sequence=' + mafft_output,
                    '-outseq=' + output,
                ],
                capture_output=True,
                text=True,
                check=True,
            )
        except subprocess.CalledProcessError as e:
            # If the bash script fails, capture the error and log the traceback
            logging.error('generate_consensus_mafft failed with error: %s', e)
            logging.error('Script stderr: %s', e.stderr)
            logging.error(traceback.format_exc())
            raise
