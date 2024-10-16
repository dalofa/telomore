"""Telomore
Script for finding and extracting telomeres from Illumina reads, which have been excluded/missed from
a de novo assembly.

Main differences from main_nanopore:
- Uses Bowtie2 to map reads
- Uses Mafft rather than lamassemble to generate consensus
- Use the quality of reads rather than number and ratio to trim ends ()
"""



# imports
import glob
import argparse
import os
import shutil
import sys

# Illumina-specific imports
from utlis.qc_reports import cons_genome_map, qc_map_illumina, finalize_log
from utlis.fasta_tools import get_chromosome, strip_fasta
from utlis.cmd_tools import map_and_sort_illumina, generate_consensus_mafft, map_and_sort
from utlis.map_tools import get_terminal_reads, get_left_soft, get_right_soft, revcomp_reads, revcomp, stich_telo, trim_by_map_illumina


def main():
    args,ref_name = get_args()
    folder_content = os.listdir()

    #
    chrom_out= ref_name + ".chrom.fa"
    get_chromosome(args.reference, chrom_out)

    # 0: Map reads and extract terminally-extending sequence
    # -----------------------------------------------------------------
    print("Mapping reads to chromsome")
    
    map_out = ref_name+".map.sort.bam"
    if not map_out in folder_content:
        if args.mode="nanopore":
            map_and_sort(chrom_out, args.fastq , map_out , args.threads)
        elif args.mode=="illumina":
            map_and_sort_illumina(chrom_out, args.fastq , map_out , args.threads)
    else:
        print("Using already identified .bam-file", map_out)

    print("Extracting terminal reads")
    # Select soft-clipped reads that extend the reference
    left_reads = "left.sam"
    right_reads = "right.sam"
    left_filt = "left_filtered"
    right_filt = "right_filtered"

    get_terminal_reads(map_out,left_reads,right_reads)
    get_left_soft(left_reads,left_filt,offset=500) 
    get_right_soft(right_reads,right_filt,offset=500)

    # 1: Generate consensus 
    # -----------------------------------------------------------------    
    
    # Generate left-consensus
    # To maintain anchor point for alignment, the reads are flipped
    # and the resulting consensus must then be flipped
    revcomp_reads("left_filtered.fastq","rev_left_filtered.fastq") # flip reads for 
    l_cons_out="rev_left_cons.fasta"
    generate_consensus_mafft("rev_left_filtered.fastq",l_cons_out)
    l_cons_final_out="left_cons.fasta"
    revcomp(l_cons_out,l_cons_final_out)


    # Generate right-consensus
    r_cons_out="right_cons.fasta"
    generate_consensus_mafft("right_filtered.fastq",r_cons_out)
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



def get_args():
    """Parse arguments"""
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
        '--mode', 
        choices=['nanopore', 'illumina'], 
        required=True, 
        help="Choose which script to run"
    )

    # Check if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # Generate a filename stripped of the .fasta/.fna/.fa extension
    ref_name = os.path.splitext(os.path.basename(args.reference))[0]
    return args, ref_name


# Run script
if __name__ == '__main__':
    main()
