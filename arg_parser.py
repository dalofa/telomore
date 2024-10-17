# arg_parser.py

import argparse
import logging
import sys

def get_args():
    """Handles parsing of arguments"""
    parser = argparse.ArgumentParser(
        description="Recover potential telomeric seq from Streptomyces Oxford Nanopore data")

    parser.add_argument(
        "-f", "--fastq", 
        type=str,
        required=True, 
        help="Path to gzipped fastq-file",
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
        "-m", '--mode', 
        choices=['nanopore', 'illumina'], 
        required=True, 
        help="Choose which mode to run"
    )

    # Check if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    
    return args

def setup_logging(log_file="telomore.log"):
    """Set-up logging"""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),  # Log file
            logging.StreamHandler(sys.stdout)  # Print to console
        ]
    )
