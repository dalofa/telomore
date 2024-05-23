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
import glob
import argparse
import os
from qc_reports import qc_map, cons_length, cons_genome_map, cons_cons_map
from fasta_tools import get_chromosome, strip_fasta
from cmd_tools import map_and_sort, train_lastDB, generate_consensus, polish
from sam_tools import get_terminal_reads, get_left_soft, get_right_soft, revcomp_reads, revcomp, stich_telo, trim_by_map
import os
import shutil

def main():
    
    # Read in arguments
    args,ref_name = get_args()
    folder_content = os.listdir()

    # Discarding nonchromosomal elements
    chrom_out= ref_name + ".chrom.fa"
    get_chromosome(args.reference, chrom_out)

    # 0: Mapping of reads and extraction of soft-clipped reads
    # -----------------------------------------------------------------
    map_out = ref_name+".map.sam"
    map_finish= map_out+".sort.bam"

    if not map_finish in folder_content:
        map_and_sort(chrom_out,args.fastq,args.threads,map_out)
        print("Mapping finished")
    else:
        print("Using already identified mapping file", map_finish)

    # Select soft-clipped reads that extend the reference
    left_reads = "left.sam"
    right_reads = "right.sam"
    left_filt = "left_filtered"
    right_filt = "right_filtered"

    get_terminal_reads(map_finish,left_reads,right_reads)
    get_left_soft(left_reads,left_filt,offset=500) 
    get_right_soft(right_reads,right_filt,offset=500)
    print("Terminals reads extracted")
    


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

    l_map_out = "left.map.sam"
    l_map_finish= l_map_out+".sort.bam"
    map_and_sort(l_map_in,"left_cons.fasta",args.threads,l_map_out)
    
    r_map_out = "right.map.sam"
    r_map_finish= r_map_out+".sort.bam"
    map_and_sort(r_map_in,"right_cons.fasta",args.threads,r_map_out)

    stitch_out = ref_name +".01.nontrimmed.cons.fasta"
    cons_log_out = ref_name + ".cons.log.txt"
    stich_telo(chrom_out,l_map_finish,r_map_finish,stitch_out,logout=cons_log_out)
    print("Consensus attached to genome")

    # 2: Trim consensus
    # ----------------------------------------------------------------- 
    trim_map = ref_name + ".01.nontrimmed.map.sam"
    trim_map_bam = trim_map + ".sort.bam"
    qc_map(stitch_out,left_reads,right_reads,trim_map,t=args.threads)
    trim_out = ref_name + ".02.trimmed.cons.fasta"
    trim_by_map(stitch_out,trim_map_bam,trim_out, cons_log=cons_log_out, cov_thres=5, ratio_thres=0.7,qual_thres=5)
    print("Consensus quality trimmed")

    # 3: QC and clean-up
    # -----------------------------------------------------------------
    qc_out = ref_name + ".02.reads.trimmed.map.sam"
    qc_map(trim_out,left_reads,right_reads,qc_out,t=args.threads)
    cons_genome_map_out = ref_name + ".02.cons.trimmed.map.sam"
    cons_genome_map("left_cons.fasta","right_cons.fasta",trim_out,cons_genome_map_out,t=args.threads)
    #cons_cons_map_out=ref_name+".03.cons_vs_cons.map"
    #cons_cons_map("left_cons.fasta","right_cons.fasta",cons_cons_map_out,t=args.threads)
    print("QC report and alignments generated")

    #rm all the files that were made
    
    # all bam-files
    MAPS_TO_REMOVE = [map_finish,l_map_finish,r_map_finish,] # trim_map_bam temporarily removed from list
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

    CONS_EXT=[".fasta",".aln"]
    for file_ext in CONS_EXT:
        os.remove(r_cons_out+file_ext)
        os.remove(l_cons_out+file_ext)
     
    os.remove(l_cons_final_out)
    os.remove("all_terminal_reads.fastq")
    #os.remove(stitch_out)

    # simplify all mapping names
    for file_name in glob.glob("*map.sam.sort*"):
        old_name = file_name
        new_name = file_name.split("map.sam.sort.")[0]+file_name.split("map.sam.sort.")[1]
        shutil.move(old_name,new_name)

    # move qc  files to QC_folder
    qc_path= ref_name + "_QC"
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
    print("Running telomore 0.1")
    main()
