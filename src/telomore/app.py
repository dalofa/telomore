"""
Telomore main application module.

Script for finding and extracting telomeres from nanopore or illumina reads,
which have been excluded from a de novo assembly.
"""

from argparse import Namespace
import logging
import os
import shutil
import traceback

from telomore._version import __version__
from telomore.utils.arg_parser import get_args, setup_logging
from telomore.utils.classes_and_small_func import Replicon
from telomore.utils.cmd_tools import (
    generate_consensus_lamassemble,
    generate_consensus_mafft,
    map_and_sort,
    map_and_sort_illumina,
    map_and_sort_illumina_cons,
    train_lastDB,
)
from telomore.utils.fasta_tools import (
    build_extended_fasta,
    extract_contig,
    get_fasta_length,
    get_linear_elements,
    strip_fasta,
)
from telomore.utils.map_tools import (
    get_left_soft,
    get_right_soft,
    get_terminal_reads,
    revcomp,
    revcomp_reads,
    stitch_telo,
    trim_by_map,
    trim_by_map_illumina,
)
from telomore.utils.qc_reports import (
    finalize_log,
    qc_map,
    qc_map_illumina,
)


def check_dependencies(required_tools: list[str] | None = None) -> None:
    """
    Check if required external dependencies are available in PATH.

    Verifies that all bioinformatics tools required by Telomore are installed
    and accessible. Logs the path to each found tool and exits with error if
    any tools are missing.

    Parameters
    ----------
    required_tools : list of str or None, optional
        List of command-line tool names to check. If None, no tools are checked.
        Common tools include: minimap2, samtools, lamassemble, mafft, bowtie2,
        lastdb, lastal, cons.

    Returns
    -------
    None
        Logs tool locations or exits if dependencies are missing.

    Raises
    ------
    SystemExit
        If any required tools are not found in PATH (exits with code 1)

    Notes
    -----
    For each tool, this function:
    - Checks if the tool is available using shutil.which()
    - Logs the full path if found
    - Collects missing tools and reports them all at once before exiting

    This ensures users know about all missing dependencies upfront rather than
    discovering them one at a time during execution.
    """
    missing_tools = []
    for tool in required_tools:
        if shutil.which(tool) is None:
            # Log missing tool
            missing_tools.append(tool)
        else:
            # Log the path to the tool
            logging.info(f'{tool}\t {shutil.which(tool)}')
    if missing_tools:
        # Log all missing tools and exit
        logging.error(f'Missing required tools: {", ".join(missing_tools)}')
        exit(1)


def entrypoint() -> None:
    """
    Entry point for the telomore command-line interface.

    Parses command-line arguments, sets up logging, and calls the main workflow.
    This function serves as the entry point defined in pyproject.toml for the
    'telomore' console script.

    Returns
    -------
    None
        Executes the main workflow or exits with error code 1 on failure.

    Raises
    ------
    SystemExit
        If argument parsing fails or an unhandled exception occurs during workflow

    Notes
    -----
    Error handling:
    - Captures all exceptions during workflow execution
    - Logs full traceback to log file
    - Exits with code 1 to signal failure to calling process

    Logging is configured before main() is called, with output to both
    console and telomore.log file (unless --quiet is specified).
    """
    args = get_args()  # Get arguments
    setup_logging(log_file='telomore.log', quiet=args.quiet)  # setup logging
    try:
        main(args)

    except Exception:
        logging.error('An error occurred during the workflow:')
        logging.error(traceback.format_exc())
        exit(1)


def main(args: Namespace) -> None:
    """
    Execute the main Telomore telomere extension workflow.

    Orchestrates the complete pipeline for extending linear contigs with
    telomeric sequences identified from unmapped reads. Processes either
    Oxford Nanopore or Illumina sequencing data based on the mode parameter.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments containing:
        - mode : str - Sequencing platform ('nanopore' or 'illumina')
        - reference : str - Path to reference genome FASTA
        - single : str - Nanopore FASTQ file (if mode='nanopore')
        - read1, read2 : str - Illumina paired FASTQ files (if mode='illumina')
        - threads : int - Number of threads for parallel operations
        - keep : bool - Whether to retain intermediate files
        - quiet : bool - Suppress console logging
        - coverage_threshold : int or None - Minimum coverage for consensus trimming
        - quality_threshold : int or None - Minimum base quality for consensus trimming

    Returns
    -------
    None
        Creates output directory with extended assemblies and QC files.

    Raises
    ------
    SystemExit
        If output folder exists, no linear contigs found, or dependencies missing

    Notes
    -----
    Workflow steps:
    1. Check external tool dependencies
    2. Identify linear contigs from reference headers
    3. Map reads to reference genome
    4. Extract terminal extending reads for each linear contig
    5. Generate consensus sequences from extending reads
    6. Align and attach consensus to contig ends
    7. Trim consensus based on read support
    8. Generate QC BAM files for manual inspection
    9. Create final assembly combining extended and unmodified contigs
    10. Clean up intermediate files (unless --keep specified)

    Platform-specific defaults:
    - Nanopore: coverage_threshold=5, quality_threshold=10
    - Illumina: coverage_threshold=1, quality_threshold=30

    Output structure: {reference_basename}_{np|ill}_telomore/
    """
    logging.info(f'Running Telomore: {__version__}')

    check_dependencies(
        [
            'minimap2',
            'samtools',
            'lamassemble',
            'mafft',
            'bowtie2',
            'lastdb',
            'lastal',
            'cons',
        ]
    )

    ref_name = os.path.splitext(os.path.basename(args.reference))[0]
    folder_content = os.listdir()

    # Create output folder
    if args.mode == 'nanopore':
        telo_folder = ref_name + '_np_telomore'
    elif args.mode == 'illumina':
        telo_folder = ref_name + '_ill_telomore'

    if os.path.isdir(telo_folder):
        logging.info('Output folder %s already exists.', telo_folder)
        exit()
    os.mkdir(telo_folder)

    # Identify linear elements
    linear_elements = get_linear_elements(args.reference)
    if not linear_elements:
        logging.info('No tagged linear elements identified')
        exit()
    logging.info('Identified the following tagged linear elements %s', linear_elements)

    # Create a list of replicon instances
    replicon_list = [Replicon(element, args.reference) for element in linear_elements]

    # 0: Map reads and extract terminally-extending sequence
    # -----------------------------------------------------------------
    logging.info('Mapping reads to assembly')

    map_out = ref_name + '_map.bam'

    # Use already existing map
    if map_out in folder_content:
        logging.info('Using already identified .bam-file %s', map_out)
    elif args.mode == 'nanopore':
        map_and_sort(
            reference=args.reference,
            fastq=args.single,
            output=map_out,
            threads=args.threads,
        )
    elif args.mode == 'illumina':
        map_and_sort_illumina(
            reference=args.reference,
            read1=args.read1,
            read2=args.read2,
            output=map_out,
            threads=args.threads,
        )

    for replicon in replicon_list:
        logging.info('\tContig %s', replicon.name)

        get_terminal_reads(
            sorted_bam_file=map_out,
            contig=replicon.name,
            loutput_handle=replicon.left_sam,
            routput_handle=replicon.right_sam,
        )
        get_left_soft(
            sam_file=replicon.left_sam, left_out=replicon.left_filt, offset=500
        )
        get_right_soft(
            sam_file=replicon.right_sam,
            contig=replicon.name,
            right_out=replicon.right_filt,
            offset=500,
        )

    # 1: Generate consensus
    # -----------------------------------------------------------------
    logging.info('Generating consensus')

    # Generate consensus
    for replicon in replicon_list:
        logging.info('\tContig %s', replicon.name)

        # GENERATE LEFT CONSENSUS
        # To maintain alignment anchor point, the reads are flipped
        # And the resulting consensus must then be flipped again
        revcomp_reads(reads_in=replicon.left_filt_fq, reads_out=replicon.revcomp_out)

        if args.mode == 'nanopore':
            db_out = ref_name + '.db'
            train_lastDB(
                args.reference, args.single, db_out, args.threads
            )  # train on entire reference
            generate_consensus_lamassemble(
                db_name=db_out, reads=replicon.revcomp_out, output=replicon.l_cons_out
            )

        elif args.mode == 'illumina':
            generate_consensus_mafft(
                reads=replicon.revcomp_out, output=replicon.l_cons_out
            )
        # flip consensus to match original orientation
        revcomp(fasta_in=replicon.l_cons_out, fasta_out=replicon.l_cons_final_out)

        # GENERATE RIGHT CONSENSUS
        # The right reads are already oriented with the anchor point
        # left-most and does therefore not need to be flipped
        if args.mode == 'nanopore':
            # A last-db should aldready exist from the left-consensus
            generate_consensus_lamassemble(
                db_name=db_out,
                reads=replicon.right_filt_fq,
                output=replicon.r_cons_final_out,
            )
        elif args.mode == 'illumina':
            generate_consensus_mafft(
                reads=replicon.right_filt_fq, output=replicon.r_cons_final_out
            )
    # 2: Extend assembly with consensus by mapping onto chromsome
    # -----------------------------------------------------------------
    logging.info('Extending assembly')

    for replicon in replicon_list:
        logging.info('\tContig %s', replicon.name)

        # Produce fasta file of just the contig to be extended
        extract_contig(
            fasta_in=replicon.org_fasta,
            contig_name=replicon.name,
            fasta_out=replicon.contig_fasta,
        )

        # Discard bases that provide alternative mapping sites
        # for the consensus to map to as Streptomyces have TIRs.
        # discard half the contig

        strip_size = int(
            get_fasta_length(
                fasta_file=replicon.contig_fasta, contig_name=replicon.name
            )
            / 2
        )
        strip_fasta(
            input_file=replicon.contig_fasta,
            output_file=replicon.trunc_left_fasta,
            x=strip_size,
            remove_from='end',
        )
        strip_fasta(
            input_file=replicon.contig_fasta,
            output_file=replicon.trunc_right_fasta,
            x=strip_size,
            remove_from='start',
        )

        if args.mode == 'nanopore':
            # Map onto the reduced reference using minimap2
            map_and_sort(
                reference=replicon.trunc_left_fasta,
                fastq=replicon.l_cons_final_out,
                output=replicon.l_map_out,
                threads=args.threads,
            )

            map_and_sort(
                reference=replicon.trunc_right_fasta,
                fastq=replicon.r_cons_final_out,
                output=replicon.r_map_out,
                threads=args.threads,
            )

        elif args.mode == 'illumina':
            # Map onto reduced reference using bowtie2
            map_and_sort_illumina_cons(
                reference=replicon.trunc_left_fasta,
                consensus_fasta=replicon.l_cons_final_out,
                output=replicon.l_map_out,
                threads=args.threads,
            )

            map_and_sort_illumina_cons(
                reference=replicon.trunc_right_fasta,
                consensus_fasta=replicon.r_cons_final_out,
                output=replicon.r_map_out,
                threads=args.threads,
            )

        # Extend the assembly using the map
        if args.mode == 'nanopore':
            cons_log_out = replicon.cons_log_np_out
        elif args.mode == 'illumina':
            cons_log_out = replicon.cons_log_ill_out

        stitch_telo(
            ref=replicon.contig_fasta,
            left_map=replicon.l_map_out,
            right_map=replicon.r_map_out,
            outfile=replicon.stitch_out,
            logout=cons_log_out,
            tmp_left=replicon.stitch_left_fasta,
            tmp_right=replicon.stitch_right_fasta,
        )

    # 3: Trim consensus using a map of terminal reads onto extended
    # assembly
    # -----------------------------------------------------------------
    logging.info('Trimming consensus based on read support')

    for replicon in replicon_list:
        # Iterate to the correct log
        if args.mode == 'nanopore':
            cons_log_out = replicon.cons_log_np_out
        elif args.mode == 'illumina':
            cons_log_out = replicon.cons_log_ill_out

        logging.info('\tContig %s', replicon.name)

        if args.mode == 'nanopore':
            # Set default values for consensus trimming if the User did not
            if args.coverage_threshold is None:
                args.coverage_threshold = 5
            if args.quality_threshold is None:
                args.quality_threshold = 10

            qc_map(
                extended_assembly=replicon.stitch_out,
                left=replicon.left_sam,
                right=replicon.right_sam,
                output_handle=replicon.trim_map,
                t=args.threads,
            )

            trim_by_map(
                untrimmed_assembly=replicon.stitch_out,
                sorted_bam_file=replicon.trim_map,
                output_handle=replicon.trim_out,
                cons_log=cons_log_out,
                cov_thres=args.coverage_threshold,
                ratio_thres=0.7,
                qual_thres=args.quality_threshold,
            )

        elif args.mode == 'illumina':
            # Set default values for consensus trimming if the User did not
            if args.coverage_threshold is None:
                args.coverage_threshold = 1
            if args.quality_threshold is None:
                args.quality_threshold = 30

            qc_map_illumina(
                extended_assembly=replicon.stitch_out,
                left_sam=replicon.left_sam,
                right_sam=replicon.right_sam,
                fastq_in1=args.read1,
                fastq_in2=args.read2,
                output_handle=replicon.trim_map,
                t=args.threads,
            )
            trim_by_map_illumina(
                untrimmed_assembly=replicon.stitch_out,
                sorted_bam_file=replicon.trim_map,
                output_handle=replicon.trim_out,
                cons_log=cons_log_out,
                cov_thres=args.coverage_threshold,
                ratio_thres=0.7,
                qual_thres=args.quality_threshold,
            )
    # 4: Generate QC files
    # -----------------------------------------------------------------
    logging.info('Generating QC map and finalizing result-log')

    for replicon in replicon_list:
        # Iterate to the correct log
        if args.mode == 'nanopore':
            cons_log_out = replicon.cons_log_np_out
        elif args.mode == 'illumina':
            cons_log_out = replicon.cons_log_ill_out

        logging.info('\tContig %s', replicon.name)

        if args.mode == 'nanopore':
            qc_map(
                extended_assembly=replicon.trim_out,
                left=replicon.left_sam,
                right=replicon.right_sam,
                output_handle=replicon.qc_out,
                t=args.threads,
            )
        if args.mode == 'illumina':
            qc_map_illumina(
                extended_assembly=replicon.trim_out,
                left_sam=replicon.left_sam,
                right_sam=replicon.right_sam,
                fastq_in1=args.read1,
                fastq_in2=args.read2,
                output_handle=replicon.qc_out,
                t=args.threads,
            )

        finalize_log(
            log=cons_log_out,
            right_fasta=replicon.stitch_right_fasta,
            left_fasta=replicon.stitch_left_fasta,
        )

    # 5: Clean-up
    # -----------------------------------------------------------------
    logging.info('Removing temporary files')

    finished_fasta = ref_name + '_telomore.fasta'
    build_extended_fasta(
        org_fasta=args.reference,
        linear_elements=linear_elements,
        replicon_list=replicon_list,
        output_handle=finished_fasta,
    )

    shutil.move(src=finished_fasta, dst=os.path.join(telo_folder, finished_fasta))

    for replicon in replicon_list:
        replicon.mv_files(telo_folder, args.mode)

    if args.keep is False:
        # rmv tmp files
        for replicon in replicon_list:
            replicon.cleanup_tmp_files()

        # rmv lastdb
        last_db_ext = ['.bck', '.des', '.par', '.prj', '.sds', '.ssp', '.suf', '.tis']

        if args.mode == 'nanopore':
            for ext in last_db_ext:
                db_file = f'{db_out}{ext}'
                os.remove(db_file)

        # remove map
        os.remove(map_out)  # map
        os.remove(f'{map_out}.bai')  # index

    logging.info('Output files moved to %s', telo_folder)
