"""Telomore
Script for finding and extracting telomeres from Illumina reads, which have been excluded/missed from
a de novo assembly.

Main differences from main_nanopore:
- Uses Bowtie2 to map reads
- Use the quality of reads rather than number and ratio to trim ends ()
"""



# imports
import glob
import argparse
import os
from qc_reports import qc_map, cons_length, cons_genome_map, cons_cons_map, qc_map_illumina, finalize_log
from fasta_tools import get_chromosome, strip_fasta
from cmd_tools import map_and_sort_illumina, train_lastDB, generate_consensus, map_and_sort
from sam_tools import get_terminal_reads, get_left_soft, get_right_soft, revcomp_reads, revcomp, stich_telo, trim_by_map_illumina
import os
import shutil
import sys

def main():
    print (__file__)
    # Read in arguments
    args,ref_name = get_args()
    folder_content = os.listdir()

    # Discarding nonchromosomal elements
    chrom_out= ref_name + ".chrom.fa"
    get_chromosome(args.reference, chrom_out)

    # 0: Mapping of reads and extraction of soft-clipped reads
    # -----------------------------------------------------------------
    map_out = ref_name+".map.sort.bam"

    if not map_out in folder_content:
        map_and_sort_illumina(chrom_out, args.fastq , map_out , args.threads)
        print("Mapping finished")
    else:
        print("Using already identified mapping file", map_out)

    # Select soft-clipped reads that extend the reference
    left_reads = "left.sam"
    right_reads = "right.sam"
    left_filt = "left_filtered"
    right_filt = "right_filtered"

    get_terminal_reads(map_out,left_reads,right_reads)
    get_left_soft(left_reads,left_filt,offset=500) 
    get_right_soft(right_reads,right_filt,offset=500)
    print("Terminal reads extracted")
    


    # 1: Generate consensus 
    # -----------------------------------------------------------------    
    db_out=ref_name+".db"
    train_lastDB(chrom_out,args.fastq,db_out,args.threads)

    # Generate left-consensus
    # To maintain anchor point for alignment, the reads are flipped
    # and the resulting consensus must then be flipped
    revcomp_reads("left_filtered.fastq","rev_left_filtered.fastq") # flip reads for 
    l_cons_out="rev_left_cons"
    generate_consensus(db_out,"rev_left_filtered.fastq",l_cons_out)
    l_cons_final_out="left_cons.fasta"
    revcomp(l_cons_out+".fasta",l_cons_final_out)


    # Generate right-consensus
    r_cons_out="right_cons"
    generate_consensus(db_out,"right_filtered.fastq",r_cons_out)
    print("Consensus generated")
    


    # Extend assembly with consensus using and mapping of the consensus
    # onto the chromosome
    # discard start or end bases that can provide alternative mapping spots for consensus
    l_map_in= ref_name + ".chrom.left.fa"
    strip_fasta(chrom_out,l_map_in,500,"end")

    r_map_in = ref_name + ".chrom.right.fa"
    
    strip_fasta(chrom_out,r_map_in ,500,"start")

    # as this is fasta and not illumina, no need to worry
    l_map_out = "left.map.sort.bam"
    map_and_sort(l_map_in,"left_cons.fasta",l_map_out,args.threads,)
    
    r_map_out = "right.map.sort.bam"
    map_and_sort(r_map_in,"right_cons.fasta",r_map_out,args.threads)

    stitch_out = ref_name +".01.nontrimmed.cons.fasta"
    cons_log_out = ref_name + ".ill.cons.log.txt"
    stich_telo(chrom_out,l_map_out,r_map_out,stitch_out,logout=cons_log_out)
    print("Consensus attached to genome")

    # 2: Trim consensus
    # ----------------------------------------------------------------- 
    print("trimming consensus")
    trim_map = ref_name + ".01.nontrimmed.map.bam"
    qc_map_illumina(stitch_out,left_reads,right_reads,trim_map,t=args.threads)
    trim_out = ref_name + ".02.trimmed.cons.fasta"
    trim_by_map_illumina(stitch_out,trim_map,trim_out, cons_log=cons_log_out, cov_thres=1, ratio_thres=0.7,qual_thres=30)
    print("Consensus quality trimmed")

    # 3: QC and clean-up
    # -----------------------------------------------------------------
    qc_out = ref_name + ".02.reads.trimmed.map.bam"
    qc_map_illumina(trim_out,left_reads,right_reads,qc_out,t=args.threads)
    cons_genome_map_out = ref_name + ".02.cons.trimmed.map.bam"
    cons_genome_map("left_cons.fasta","right_cons.fasta",trim_out,cons_genome_map_out,t=args.threads)
    #cons_cons_map_out=ref_name+".03.cons_vs_cons.map"
    #cons_cons_map("left_cons.fasta","right_cons.fasta",cons_cons_map_out,t=args.threads)
    finalize_log(cons_log_out,"tmp.right.fasta","tmp.left.fasta")
    print("QC report and alignments generated")

    #rm all the files that were made
    if args.keep==False:
    # all bam-files
        MAPS_TO_REMOVE = [map_out,l_map_out,r_map_out,] # trim_map_bam temporarily removed from list
        for file in MAPS_TO_REMOVE:
            os.remove(file)
            os.remove(file+".bai")

        DATABASE_EXT=[".bck",".des",".par",".prj",".sds",".ssp",".suf",".tis"]
        for file_ext in DATABASE_EXT:
            os.remove(db_out+file_ext)

        TERMINAL_MAP_ENDINGS=[".sam",".fastq"]
        for file_ext in TERMINAL_MAP_ENDINGS:
            os.remove(left_filt+file_ext)
            os.remove(right_filt+file_ext)
            if file_ext==".fastq":
                os.remove("rev_"+left_filt+file_ext)
        
        os.remove(left_reads)
        os.remove(right_reads)

        # Check if file exists to avoid crashing when no alignment was used to produce
        # alignment
        CONS_EXT=[".fasta",".aln"]
        for file_ext in CONS_EXT:
            if os.path.isfile(r_cons_out+file_ext):
                os.remove(r_cons_out+file_ext)
            if os.path.isfile(r_cons_out+file_ext):
                os.remove(l_cons_out+file_ext)
        
        #remove temp files not needed in the end
        os.remove(l_cons_final_out)
        os.remove(l_map_in)
        os.remove(r_map_in)
        os.remove("tmp.left.fasta")
        os.remove("tmp.right.fasta")
        os.remove("all_terminal_reads.fastq")
        #os.remove(stitch_out)

    # simplify all mapping names
    for file_name in glob.glob("*map.sam.sort*"):
        old_name = file_name
        new_name = file_name.split("map.sam.sort.")[0]+file_name.split("map.sam.sort.")[1]
        shutil.move(old_name,new_name)

    # move qc  files to QC_folder
    qc_path= ref_name + "_ill" + "_QC"
    os.mkdir(qc_path)

    for QC_FILE in glob.glob(ref_name+".0*"):
        
        if "02.trimmed.cons.fasta" in QC_FILE:
            shutil.copyfile(QC_FILE,os.path.join(qc_path,QC_FILE))
        else:
            shutil.move(QC_FILE, os.path.join(qc_path,QC_FILE))


    # outcommented and replaced by the glob.glob utilizing code above
    # QC_FILES = [qc_out,cons_genome_map_out,cons_cons_map_out]
    # for file in QC_FILES:
    #     mv1 = file + ".sort.bam"
    #     mv2 = mv1 + ".bai"
    #     shutil.move(mv1, os.path.join(qc_path,mv1))
    #     shutil.move(mv2, os.path.join(qc_path,mv2))
    # shutil.move(cons_log_out,os.path.join(qc_path,cons_log_out))


def get_args():
    """Parses arguments"""
    parser = argparse.ArgumentParser(
        description="Recover potential telomeric seq from Streptomyces Oxford Nanopore data")

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

    # Check if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # Validate input options
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
    ref_name = args.reference.split(".")[0]
    if "\\" in ref_name:
        ref_name = ref_name.split("\\")[-1]
    elif "/" in ref_name:
        ref_name = ref_name.split("/")[-1]

    return args, ref_name


# Run script
if __name__ == '__main__':
    print("Running telomore 0.1 in Illumina mode")
    main()
