# arg_parser.py

import argparse
import logging
import sys

def get_args():
    """Handles parsing of arguments"""
    parser = argparse.ArgumentParser(
        description= """Recover potential telomeric sequences from Streptomyces genomes using Oxford Nanopore or illumina sequence data.
- INPUT: Takes a fasta file with linear contigs tagged "linear" in the fasta_header and sequencing reads in fq.gz format.
- OUTPUT: An extended assembly is written to basename.02.trimmed.fasta and QC-maps are written to a folder 
  named basename_seqtype_QC.
- LOG: A run-log is written to telomore.log and a result-log is written to basename.seqtype.cons.log.txt.
""",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "-m", '--mode', 
        choices=['nanopore', 'illumina'], 
        required=True, 
        help="""Choose which mode to run.
        --mode=nanopore takes a single read-file, specified using --single
        --mode=illumina-mode takes two read-files with one mate-pair in each, specified using --read1 and --read2"""
    )
    parser.add_argument(
        "--single", 
        type=str,
        help="Path to a single gzipped nanopore fastq-file",
    )
    parser.add_argument(
        "--read1", 
        type=str,
        help="Path to gzipped illumina mate1 fastq-file",
    )
    parser.add_argument(
        "--read2", 
        type=str,
        help="Path to gzipped illumina mate2 fastq-file",
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
