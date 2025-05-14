"""
Argument parser for the telomore tool.
"""

import argparse
from argparse import Namespace
import logging
import sys

def get_args() -> Namespace:
    """Handles parsing of arguments"""
    parser = argparse.ArgumentParser(
        description= """Telomore: A tool to recover potential telomeric sequences from Streptomyces genomes.

This tool processes sequencing data from Oxford Nanopore or Illumina platforms to extend assemblies and generate QC reports.

INPUT:
- For Nanopore mode (--mode=nanopore): Provide a single gzipped FASTQ file using --single.
- For Illumina mode (--mode=illumina): Provide two gzipped FASTQ files using --read1 and --read2.
- A reference genome file in FASTA format is required for both modes (--reference).

OUTPUT:
- Extended assembly written to basename.02.trimmed.fasta (basename is the name of the input file without the extension).
- QC reports saved in a folder named basename_seqtype_QC.
- Logs are written to telomore.log and basename.seqtype.cons.log.txt. in basename_seqtyope_QC.

OPTIONS:
- Specify the number of threads to use with --threads (default: 1).
- Use --keep to retain intermediate files (default: False).
- Use --quiet to suppress console logging.

EXAMPLES:
1. Nanopore mode:
   python arg_parser.py --mode=nanopore --single reads.fastq.gz --reference genome.fasta

2. Illumina mode:
   python arg_parser.py --mode=illumina --read1 read1.fastq.gz --read2 read2.fastq.gz --reference genome.fasta
""",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "-m", '--mode', 
        choices=['nanopore', 'illumina'], 
        required=True, 
        help="""Choose which mode to run.
        --mode=nanopore takes a single read-file, specified using --single
        --mode=illumina-mode takes two read-files specified using --read1 and --read2"""
    )
    parser.add_argument(
        "--single", 
        type=str,
        help="Path to a single gzipped nanopore fastq-file",
    )
    parser.add_argument(
        "--read1", 
        type=str,
        help="Path to gzipped illumina read1 fastq-file",
    )
    parser.add_argument(
        "--read2", 
        type=str,
        help="Path to gzipped illumina read2 fastq-file",
    )
    parser.add_argument(
        "-r", "--reference", 
        type=str, 
        required=True, 
        help="Path to reference file (.fasta, .fna, or .fa)"
    )
    parser.add_argument(
        "-t", "--threads", 
        type=int, 
        default=1,
        help="Threads to use. Default is 1"
    )
    parser.add_argument(
        "-k", "--keep", 
        action='store_true', 
        help="Flag to keep intermediate files. Default is False"
    )
    parser.add_argument(
        "-q", "--quiet", 
        action='store_true', 
        help="Set logging to quiet."
    ) 

    # Check if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    
    if args.mode == 'illumina':
        if not (args.read1 and args.read2):
            parser.error('Illumina mode requires two FASTQ files, specified by --read1 and --read2')
    elif args.mode == 'nanopore':
        if not args.single:
            parser.error('Nanopore mode takes one collected FASTQ file, specified by --single')
    return args


def setup_logging(log_file="telomore.log",quiet:bool=False):
    """Set-up logging"""

    if quiet==True:
        handlers_to_use =[
            logging.FileHandler(log_file),  # Log file
        ]
    else:
        handlers_to_use =[
            logging.FileHandler(log_file),  # Log file
            logging.StreamHandler(sys.stdout)  # Print to console
        ]
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(message)s",
        handlers=handlers_to_use
    )
