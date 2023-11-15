"""Telomore
Script for finding and extracting telomeres from Nanopore reads, which have been excluded/missed from
a de novo assembly.
"""

# To implement
#0: Wrong joining of consensus to original genome
#1. QC reports: cons lenght, qc reads, qc cons
#2. Trimming according to coverage
#3. Filtering of reads for base quality
#4. Check for existing mapping file
#5. How to provide arguments
#6. Clean-up
#X. after telomore: clustering of telomers

# imports
import argparse
import os
from qc_reports import qc_map, cons_length, cons_genome_map, cons_cons_map
from fasta_tools import get_chromosome, attach_seq
from cmd_tools import map_and_sort, train_lastDB, generate_consensus, polish
from sam_filters import get_terminal_reads, get_left_soft, get_right_soft, revcomp_reads, revcomp
import shutil
import os


def main():
    print("Running telomore 0.1")
    print("Working in",os.getcwd())
    
    # prepare STUFF
    args,ref_name = get_args()
    folder_cont = os.listdir()

    # 1. Sort away non-chromsomal contigs
    chrom_out = ref_name + "_chrom.fa"
    get_chromosome(args.reference, chrom_out)

    # 2. Map reads
    map_out = ref_name+".map.sam"
    map_finish= map_out+".sort.bam"
    if not map_finish in folder_cont:
        map_and_sort(chrom_out,args.fastq,args.threads,map_out)
        print("Mapping finished")
    else:
        print("Using already identified mapping file", map_out+".sort.bam")

    # 3. Filtering of reads
    get_terminal_reads(map_finish,"left.sam","right.sam")

    get_left_soft("left.sam","left_filtered",offset=100)
    
    get_right_soft("right.sam","right_filtered",offset=100)
    
    print("Terminal reads filtered and extracted")

    # 4. Generate consensus
    #       4.1 Train lastDB
    db_out=ref_name+"_db"
    train_lastDB(chrom_out,args.fastq,db_out,args.threads)
    # this needs to happen for both sides min vennebror i kristihimmelfart
    
    #       4.2 Produces left-consensus (which requires flipping reads and then flipping the resulting consensus)
    revcomp_reads("left_filtered.fastq","rev_left_filtered.fastq") # flip reads for 
    generate_consensus(db_out,"rev_left_filtered.fastq","rev_left_cons")
    revcomp("rev_left_cons.fasta","left_cons.fasta")

    #       4.3 Produce right-consensus
    generate_consensus(db_out,"right_filtered.fastq","right_cons")

    # 5. Attach and polish
    ref_att_cons=ref_name+"_att_cons.fa"

    attach_seq("left_cons.fasta","right_cons.fasta",chrom_out,ref_att_cons,offset=100)
    
    polish(ref_att_cons,args.fastq,"medaka",threads=args.threads)
    ref_att_cons_polish=ref_name+"cons_polished.fa"
    #move file baaaack, this aint working, shoudl be now
    shutil.copyfile("medaka/consensus.fasta", ref_att_cons_polish)

    # 6. Evaluate consensus
    # QC report:
    
    
    qc_map(ref_att_cons_polish,"left.sam","right.sam",ref_name+".qc.map.sam",t=args.threads)
    cons_genome_map("left_cons.fasta","right_cons.fasta",ref_att_cons_polish,ref_name+".cons.map",t=args.threads)
    cons_cons_map("left_cons.fasta","right_cons.fasta",ref_name+".cons_vs_cons.map.sam")
    cons_length("all_cons.fasta","cons_lenghts.log.txt")
    

    # Bam-file of org. cons and reads against full polished genome


    #get reads and assemble reference

    
    # parse args
    # get chromosome
    # map reads
    # filter reads
    # extract overhanging reads


    # identify overhangs at termini
    # create consensus
    # evaluate consensus: reads supporting it and its quality at each position
    # polish consensus
    # evaluate reads supporting consensus and its quality


# parse arguments
def get_args():
    """Parses arguments"""
    parser = argparse.ArgumentParser(
    description="Recover ssDNA from Streptomyces Oxford Nanopore data")

    parser.add_argument("fastq",type=str,help="path to gzipped fastq-file")
    parser.add_argument("reference",type=str,help="path to gzipped fastq-file")
    parser.add_argument("threads",type=int,help="Threads to use")

    args=parser.parse_args()

    # generate a filename stripped of the .fasta
    ref_name = args.reference.split(".")[-2]
    if "\\" in ref_name:
        ref_name = ref_name.split("\\")[-1]
    elif "/" in ref_name:
        ref_name = ref_name.split("/")[-1]

    return args,ref_name


# Run script
if __name__ == '__main__':
    main()