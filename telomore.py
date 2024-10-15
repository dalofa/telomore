"""Telomore
Script for finding and extracting telomeres from Nanopore OR Illumina reads, which have been excluded/missed from
a de novo assembly.

This serves as a wrapper for either the illumina or nanopore workflows.
"""

import argparse
import subprocess
import os
import sys

def parse_args():
    """Parses arguments for the wrapper script."""
    parser = argparse.ArgumentParser(description="Telomore wrapper for Nanopore and illumina mode")

    # Common arguments from both scripts
    parser.add_argument(
        "-f", "--fastq", 
        type=str, 
        help="Path to gzipped fastq-file"
    )
    parser.add_argument(
        "-r", "--reference", 
        type=str, 
        help="Path to reference file (.fasta, .fna, or .fa)"
    )
    parser.add_argument(
        "-d", "--directory", 
        type=str, 
        help="Path to directory containing both the reference (.fasta, .fna, or .fa) and reads (.fq.gz or .fastq.gz)"
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

    # Mode argument to decide which script to run
    parser.add_argument(
        '--mode', 
        choices=['nanopore', 'illumina'], 
        required=True, 
        help="Choose which script to run"
    )

    # Check if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()

def validate_and_prepare_args(args):
    """Validate and prepare arguments, especially when using the --directory option."""
    
    if args.directory:
        # If directory is provided, locate the reference and fastq files
        ref_file = None
        fastq_file = None
        for file in os.listdir(args.directory):
            if file.endswith((".fasta", ".fna", ".fa")):
                ref_file = os.path.join(args.directory, file)
            elif file.endswith(".fq.gz") or file.endswith(".fastq.gz"):
                fastq_file = os.path.join(args.directory, file)

        if not ref_file or not fastq_file:
            sys.stderr.write("Error: Could not find both .fasta/.fna/.fa (reference) and .fq.gz/.fastq.gz (reads) files in the directory.\n")
            sys.exit(1)

        args.reference = ref_file
        args.fastq = fastq_file

    elif not args.fastq or not args.reference:
        sys.stderr.write("Error: Either provide both --fastq and --reference or use --directory to specify a folder containing them.\n")
        sys.exit(1)

    # Generate a filename stripped of the .fasta/.fna/.fa extension
    ref_name = args.reference.split(".")[-2]
    if "\\" in ref_name:
        ref_name = ref_name.split("\\")[-1]
    elif "/" in ref_name:
        ref_name = ref_name.split("/")[-1]

    return ref_name

def main():
    args = parse_args()

    # Validate arguments and prepare necessary values
    ref_name = validate_and_prepare_args(args)

    # Build the argument list to pass to the respective script
    script_args = [
        '--fastq', args.fastq,
        '--reference', args.reference,
        '--threads', str(args.threads)
    ]
    
    if args.keep:
        script_args.append('--keep')

    # Based on the mode argument, call the respective script
    basedir = os.path.dirname(__file__) # nessesary to find the location of the python scripts
    if args.mode == 'nanopore':
        subprocess.run([sys.executable, os.path.join(basedir,'main_nanopore.py')] + script_args)
    elif args.mode == 'illumina':
        subprocess.run([sys.executable, os.path.join(basedir,'main_illumina.py')] + script_args)

if __name__ == "__main__":
    main()