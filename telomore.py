"""Telomore
Script for finding and extracting telomeres from nanopore or illumina reads, which have been excluded from
a de novo assembly.
"""
# imports
import json
import glob
import os
import shutil
import logging
import traceback
from utils.arg_parser import get_args, setup_logging
from utils.fasta_tools import strip_fasta,get_fasta_length,get_linear_elements, extract_contig
from utils.map_tools import get_terminal_reads, get_left_soft, get_right_soft, revcomp_reads, revcomp, stich_telo, get_contig_map
from utils.qc_reports import finalize_log

# Nanopore-sepcific imports
from utils.cmd_tools import map_and_sort, generate_consensus_lamassemble, train_lastDB
from utils.qc_reports import qc_map
from utils.map_tools import trim_by_map

# Illumina-specific imports
from utils.qc_reports import qc_map_illumina
from utils.cmd_tools import map_and_sort_illumina, generate_consensus_mafft
from utils.map_tools import trim_by_map_illumina


def main():
    
    args = get_args()

    # Generate a filename stripped of the .fasta/.fna/.fa extension
    ref_name = os.path.splitext(os.path.basename(args.reference))[0]
    folder_content = os.listdir()

    logging.info(f"Running Telomore: 0.3 in {args.mode} mode")

    # Identify linear elements
    linear_elements = get_linear_elements(args.reference)
    
    # check if no linear element
    if not linear_elements:
         logging.info("No tagged linear elements identified")
         exit()

    # Intialize dict for keeping track of files associated with each replicon
    replicon_dict = {replicon: {"org_file": args.reference,} for replicon in linear_elements}

    logging.info(f"Identified the following tagged linear elements {linear_elements}")
    # 0: Map reads and extract terminally-extending sequence
    # -----------------------------------------------------------------

    logging.info("Mapping reads to assembly")

    map_out = ref_name+"_map.bam"

    if not map_out in folder_content:

        if args.mode=="nanopore":
            map_and_sort(args.reference, args.single , map_out , args.threads)
        elif args.mode=="illumina":
            map_and_sort_illumina(args.reference, args.read1, args.read2, map_out , args.threads)
    else:
        logging.info(f"Using already identified .bam-file {map_out}")

    # get map for each contig
    for replicon in linear_elements:
        contig_map_out = replicon + "_map.bam"

        get_contig_map(bam_in=map_out,
                       contig_name=replicon,
                       output_path = contig_map_out) 
        # save files produced to dict
        replicon_dict[replicon]["map"]=map_out
        replicon_dict[replicon]["contig_map"]=contig_map_out
    logging.info("Extracting terminal reads")

    # Extract extending terminal reads for each contig via the contig-specific map
    
    for replicon in replicon_dict.keys():
         logging.info(f"\tContig {replicon}")
         map = replicon_dict[replicon]["contig_map"] # should be contig_map?

         left_sam = replicon + "_left.sam"
         right_sam = replicon + "_right.sam"
         left_filt = replicon + "_left_filtered"
         right_filt = replicon + "_right_filtered"

         get_terminal_reads(map,left_sam,right_sam)
         get_left_soft(sam_file = left_sam,
                       left_out = left_filt,
                       offset=500)
         get_right_soft(sam_file = right_sam,
                       right_out = right_filt,
                       offset=500)
         
         # add files to dict to keep track of them
         replicon_dict[replicon]["left_sam"]=left_sam
         replicon_dict[replicon]["right_sam"]=right_sam
         replicon_dict[replicon]["right_filt_sam"]=right_filt + ".sam"
         replicon_dict[replicon]["left_filt_sam"]=left_filt +".sam"
         replicon_dict[replicon]["right_filt_fq"]=right_filt + ".fastq"
         replicon_dict[replicon]["left_filt_fq"]=left_filt +".fastq"

    # 1: Generate consensus
    # -----------------------------------------------------------------    
    logging.info("Generating consensus")

    # Generate left-consensus
    for replicon in replicon_dict.keys():
        logging.info(f"\tContig {replicon}")
        # To maintain alignment anchor point, the reads are flipped
        # And the resulting consensus must then be flipped again
        left_filt_reads = replicon_dict[replicon]["left_filt_fq"]
        revcomp_out = "rev_" + left_filt_reads
        revcomp_reads(left_filt_reads, revcomp_out)

        l_cons_out ="rev_" + replicon + "_left_cons.fasta"
        l_cons_final_out = replicon + "_left_cons.fasta"
        
        if args.mode=="nanopore":
                db_out=ref_name+".db"
                train_lastDB(args.reference,args.single,db_out,args.threads) # train on entire reference
                generate_consensus_lamassemble(db_out,revcomp_out,l_cons_out)
                replicon_dict[replicon]["last_db"]=db_out

        elif args.mode=="illumina":
                generate_consensus_mafft(revcomp_out,l_cons_out)

        # flip consensus to match original orientation
        revcomp(l_cons_out,l_cons_final_out)

        # add files to dict
        replicon_dict[replicon]["revcomp_out"]=revcomp_out
        replicon_dict[replicon]["l_cons_out"]=l_cons_out
        replicon_dict[replicon]["l_cons_final_out"]=l_cons_final_out

    for replicon in replicon_dict.keys():
        
        # The right reads are already oriented with the anchor point
        # left-most and does therefore not need to be flipped
        right_filt_reads = replicon_dict[replicon]["right_filt_fq"]
        r_cons_final_out = replicon + "_right_cons.fasta"
        
        if args.mode=="nanopore":
                # A last-db should aldready exist from the left-consensus
                generate_consensus_lamassemble(db_out,right_filt_reads,r_cons_final_out)
        elif args.mode=="illumina":
                generate_consensus_mafft(right_filt_reads,r_cons_final_out)

        # add files to dict
        replicon_dict[replicon]["r_cons_final_out"]=r_cons_final_out
    
    # 2: Extend assembly with consensus by mapping onto chromsome
    # -----------------------------------------------------------------    
    logging.info("Extending assembly")
    
    for replicon in replicon_dict.keys():
        logging.info(f"\tContig {replicon}")
        # Produce fasta file of just the contig to be extended
        org_fasta = replicon_dict[replicon]["org_file"]
        rep_left_cons = replicon_dict[replicon]["l_cons_final_out"]
        rep_right_cons = replicon_dict[replicon]["r_cons_final_out"]

        tmp_cons_left= replicon + "_left.fasta"
        tmp_cons_right= replicon + "_right.fasta"
        contig_fasta_out = replicon + ".fasta"
        extract_contig(fasta_in = org_fasta,
                       contig_name=replicon,
                       fasta_out=contig_fasta_out)
        # Discard bases that provide alternative mapping sites
        # For the consensus to map to as Streptomyces have TIRs.
        l_map_in = replicon + "_left.fa"
        r_map_in = replicon + "_right.fa"
        # discard half the contig
        strip_size=int(get_fasta_length(contig_fasta_out,replicon)/2)
        strip_fasta(contig_fasta_out,l_map_in,strip_size,"end")
        strip_fasta(contig_fasta_out,r_map_in ,strip_size,"start")
        
        # Map onto reduced reference
        l_map_out = replicon + "_left_map.bam"
        map_and_sort(reference = l_map_in,
                     fastq = rep_left_cons,
                     output = l_map_out,
                     threads = args.threads)
        
        r_map_out = replicon + "_right_map.bam"
        map_and_sort(reference = r_map_in,
                     fastq = rep_right_cons,
                     output = r_map_out,
                     threads = args.threads)
        
        # Extend the assembly using the map
        stitch_out = replicon +"_telomore_untrimmed.fasta"
        
        if args.mode=="nanopore":
                cons_log_out= replicon + "_telomore_ext_np.log"
        elif args.mode=="illumina":
            cons_log_out = replicon + "_telomore_ill_ext.log"

        stich_telo(ref = contig_fasta_out,
                   left_map = l_map_out,
                   right_map = r_map_out,
                   outfile = stitch_out,
                   logout=cons_log_out,
                   tmp_left=tmp_cons_left,
                   tmp_right=tmp_cons_right)

        # Add files to dict
        replicon_dict[replicon]["l_map_in"]=l_map_in
        replicon_dict[replicon]["r_map_in"]=r_map_in
        replicon_dict[replicon]["stitch_out"]=stitch_out
        replicon_dict[replicon]["cons_log_out"]=cons_log_out
        replicon_dict[replicon]["contig_fasta_out"]=contig_fasta_out
        replicon_dict[replicon]["l_map_out"]=l_map_out
        replicon_dict[replicon]["r_map_out"]=r_map_out
        replicon_dict[replicon]["tmp_cons_left"]=tmp_cons_left
        replicon_dict[replicon]["tmp_cons_right"]=tmp_cons_right
        
    # 3: Trim consensus using a map of terminal reads onto extended
    # assembly
    # ----------------------------------------------------------------- 
    logging.info("Trimming consensus based on read support")
    


    for replicon in replicon_dict.keys():
        logging.info(f"\tContig {replicon}")
        untrimmed_fasta = replicon_dict[replicon]["stitch_out"]
        cons_log_out = replicon_dict[replicon]["cons_log_out"]
        
        org_left_map = replicon_dict[replicon]["left_sam"]
        org_right_map = replicon_dict[replicon]["right_sam"]
        

        trim_map = replicon + "_telomore_untrimmed.bam"
        trim_out = replicon + "_telomore_extended.fasta"


        if args.mode=="nanopore":
            qc_map(extended_assembly=untrimmed_fasta,
                   left = org_left_map,
                   right= org_right_map,
                   output_handle=trim_map,
                   t=args.threads)  
            trim_by_map(untrimmed_assembly=untrimmed_fasta,
                        sorted_bam_file=trim_map,
                        output_handle=trim_out,
                        cons_log=cons_log_out,
                        ratio_thres=0.7,
                        qual_thres=5)
              
        elif args.mode=="illumina":
            qc_map_illumina(extended_assembly=untrimmed_fasta,
                             left=org_left_map,
                             right=org_right_map,
                             fastq_in1=args.read1,
                             fastq_in2=args.read2,
                             output_handle=trim_map,
                             t=args.threads)
            trim_by_map_illumina(untrimmed_assembly=untrimmed_fasta,
                                 sorted_bam_file = trim_map,
                                 cons_log = cons_log_out,
                                 cov_thres=1,
                                 ratio_thres=0.7,
                                 qual_threshold=30)
        # Add file to dict
        replicon_dict[replicon]["trim_map"]=trim_map
        replicon_dict[replicon]["final_assembly"]=trim_out

    
    # 4: Generate QC files
    # -----------------------------------------------------------------
    logging.info("Generating QC map and finalizing result-log")

    for replicon in replicon_dict.keys():
        logging.info(f"\tContig {replicon}")
        qc_out = replicon + "_telomore_QC.bam"
        tmp_cons_left = replicon_dict[replicon]["tmp_cons_left"]
        tmp_cons_right = replicon_dict[replicon]["tmp_cons_left"]
        
        final_assembly = replicon_dict[replicon]["final_assembly"]
        cons_log_out = replicon_dict[replicon]["cons_log_out"]

        if args.mode=="nanopore":
            qc_map(extended_assembly=final_assembly,
                            left=org_left_map,
                            right=org_right_map,
                            output_handle=qc_out,
                            t=args.threads)
        if args.mode=="illumina":
            qc_map_illumina(extended_assembly=final_assembly,
                             left=org_left_map,
                             right=org_right_map,
                             fastq_in1=args.read1,
                             fastq_in2=args.read2,
                             output_handle=qc_out,
                             t=args.threads)
            
        finalize_log(log = cons_log_out,
                     right_fasta = tmp_cons_left,
                     left_fasta = tmp_cons_right)
        # add to dict
        replicon_dict[replicon]["qc_out"]=qc_out
        
    # 5: Clean-up
    # -----------------------------------------------------------------
    logging.info("Removing temporary files")
    
    # move qc  files to QC_folder
    if args.mode=="nanopore":
        telo_folder=ref_name + "_np_telomore_0"
    elif args.mode=="illumina":
        telo_folder=ref_name + "_ill_telomore_0"

    # make sure you dont have to move the folder each time
    counter=0
    while os.path.isdir(telo_folder):
        telo_folder = "_".join(telo_folder.split("_")[:-1])
        telo_folder = telo_folder + "_" + str(counter)
        counter+=1
    
    os.mkdir(telo_folder)

    for replicon, replicon_feat in replicon_dict.items():
        # Note the non-tmp files
        # Delete everything else
        keys_to_save = ["final_assembly",
        "cons_log_out",
        "qc_out"]

        input_files = ["org_file"]

        for key in replicon_feat.keys():
            if key in keys_to_save:
                file_path = replicon_dict[replicon][key]
                shutil.move(src = file_path,
                            dst = os.path.join(telo_folder,file_path))
                if key=="qc_out":
                    map_index = replicon_dict[replicon][key] +".bai"
                    shutil.move(src = map_index,
                            dst = os.path.join(telo_folder,map_index))
            elif key in input_files:
                continue
            elif "map" in key and args.keep==False: #rm map-indices
                file_path = replicon_dict[replicon][key]
                map_index = file_path + ".bai"
                if os.path.isfile(file_path):
                     os.remove(file_path)
                if os.path.isfile(map_index):
                     os.remove(map_index)
            elif "cons" in key and args.keep==False: #rm map-indices
                file_path = replicon_dict[replicon][key]
                alignment = file_path + ".aln"
                if os.path.isfile(file_path):
                     os.remove(file_path)
                if os.path.isfile( alignment):
                     os.remove(alignment)
            elif key=="last_db" and args.keep==False: #rm all LAST dbs
                DATABASE_EXT=[".bck",".des",".par",".prj",".sds",".ssp",".suf",".tis"]
                for extension in DATABASE_EXT:
                     file_path = replicon_dict[replicon][key]+extension
                     if os.path.isfile(file_path):
                        os.remove(file_path)
            elif args.keep==False:
                file_path = replicon_dict[replicon][key]
                if os.path.isfile(file_path):
                    os.remove(file_path)

    logging.info(f"Output files moved to {telo_folder}")

# Run script
if __name__ == '__main__':
    setup_logging()

    try:
        main()
    
    except Exception as e:
        logging.error("An error occurred during the workflow:")
        logging.error(traceback.format_exc()) 
