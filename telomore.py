"""
Telomore: Script for finding and extracting telomeres from nanopore or illumina reads, which have been 
excluded from a de novo assembly.
"""

# imports
import os
import shutil
import logging
import traceback
from utils import *

def main(args):
    # Generate a filename stripped of the .fasta/.fna/.fa extension
    ref_name = os.path.splitext(os.path.basename(args.reference))[0]
    folder_content = os.listdir()

    logging.info(f"Running Telomore: 0.4 [{get_git_commit_hash()}] in {args.mode} mode")
    
    # Create output folder
    if args.mode=="nanopore":
        telo_folder=ref_name + "_np_telomore"
    elif args.mode=="illumina":
        telo_folder=ref_name + "_ill_telomore"

    if os.path.isdir(telo_folder):
        logging.info(f"Output folder {telo_folder} already exsists.")
        exit()
    os.mkdir(telo_folder)

    # Identify linear elements
    linear_elements = get_linear_elements(args.reference)
    if not linear_elements:
         logging.info("No tagged linear elements identified")
         exit()
    logging.info(f"Identified the following tagged linear elements {linear_elements}")
    # Intialize dict for keeping track of files associated with each replicon
    # Assuming linear_elements is defined as a list of strings

    # Create a list of replicon instances
    replicon_list = [Replicon(element, args.reference) for element in linear_elements]   

    # 0: Map reads and extract terminally-extending sequence
    # -----------------------------------------------------------------
    logging.info("Mapping reads to assembly")

    map_out = ref_name+"_map.bam"

    # Use already existing map
    if map_out in folder_content:
         logging.info(f"Using already identified .bam-file {map_out}")
    elif args.mode=="nanopore":
         map_and_sort(reference= args.reference, 
                      fastq = args.single,
                      output = map_out,
                      threads = args.threads)
    elif args.mode=="illumina":
         map_and_sort_illumina(reference = args.reference,
                               read1 = args.read1,
                               read2 = args.read2,
                               output=map_out,
                               threads = args.threads)
    
    for replicon in replicon_list:

         logging.info(f"\tContig {replicon.name}")

         get_terminal_reads(sorted_bam_file=map_out,
                            contig=replicon.name,
                            loutput_handle=replicon.left_sam,
                            routput_handle=replicon.right_sam)
         get_left_soft(sam_file = replicon.left_sam,
                       left_out = replicon.left_filt,
                       offset=500)
         get_right_soft(sam_file = replicon.right_sam,
                        contig=replicon.name,
                        right_out = replicon.right_filt,
                        offset=500)
         
    # 1: Generate consensus
    # -----------------------------------------------------------------    
    logging.info("Generating consensus")

    # Generate consensus
    for replicon in replicon_list:
        logging.info(f"\tContig {replicon.name}")

        # GENERATE LEFT CONSENSUS
        # To maintain alignment anchor point, the reads are flipped
        # And the resulting consensus must then be flipped again
        revcomp_reads(reads_in = replicon.left_filt_fq,
                      reads_out = replicon.revcomp_out)
        
        if args.mode=="nanopore":
                db_out=ref_name+".db"
                train_lastDB(args.reference,
                             args.single,
                             db_out,
                             args.threads) # train on entire reference
                generate_consensus_lamassemble(db_name = db_out,
                                               reads = replicon.revcomp_out,
                                               output = replicon.l_cons_out)

        elif args.mode=="illumina":
                generate_consensus_mafft(reads = replicon.revcomp_out,
                                         output = replicon.l_cons_out)
        # flip consensus to match original orientation
        revcomp(fasta_in = replicon.l_cons_out,
                fasta_out = replicon.l_cons_final_out)

        # GENERATE RIGHT CONSENSUS
        # The right reads are already oriented with the anchor point
        # left-most and does therefore not need to be flipped
        if args.mode=="nanopore":
                # A last-db should aldready exist from the left-consensus
                generate_consensus_lamassemble(db_name = db_out,
                                               reads = replicon.right_filt_fq,
                                               output = replicon.r_cons_final_out)
        elif args.mode=="illumina":
                generate_consensus_mafft(reads = replicon.right_filt_fq,
                                         output = replicon.r_cons_final_out)
    # 2: Extend assembly with consensus by mapping onto chromsome
    # -----------------------------------------------------------------    
    logging.info("Extending assembly")
    
    for replicon in replicon_list:
        logging.info(f"\tContig {replicon.name}")
        
        # Produce fasta file of just the contig to be extended
        # rep_left_cons = replicon_dict[replicon]["l_cons_final_out"]
        # rep_right_cons = replicon_dict[replicon]["r_cons_final_out"]

        # tmp_cons_left= replicon + "_left.fasta"
        # tmp_cons_right= replicon + "_right.fasta"
        extract_contig(fasta_in = replicon.org_fasta,
                       contig_name = replicon.name,
                       fasta_out=replicon.contig_fasta)
        
        # Discard bases that provide alternative mapping sites
        # For the consensus to map to as Streptomyces have TIRs.
        # discard half the contig
        
        strip_size=int(get_fasta_length(fasta_file = replicon.contig_fasta,
                                        contig_name = replicon.name)/2)
        strip_fasta(input_file = replicon.contig_fasta,
                    output_file=replicon.trunc_left_fasta,
                    x = strip_size,
                    remove_from= "end")
        strip_fasta(input_file = replicon.contig_fasta,
                    output_file = replicon.trunc_right_fasta,
                    x = strip_size,
                    remove_from="start")
        if args.mode=="nanopore":
            # Map onto reduced reference
            map_and_sort(reference = replicon.trunc_left_fasta,
                        fastq = replicon.l_cons_final_out,
                        output = replicon.l_map_out,
                        threads = args.threads)
            
            map_and_sort(reference = replicon.trunc_right_fasta,
                        fastq = replicon.r_cons_final_out,
                        output = replicon.r_map_out,
                        threads = args.threads)
            
        elif args.mode=="illumina":
            # Map onto reduced reference using bowtie2
            map_and_sort_illumina_cons(reference = replicon.trunc_left_fasta,
                        consensus_fasta= replicon.l_cons_final_out,
                        output = replicon.l_map_out,
                        threads = args.threads)
            
            map_and_sort_illumina_cons(reference = replicon.trunc_right_fasta,
                        consensus_fasta = replicon.r_cons_final_out,
                        output = replicon.r_map_out,
                        threads = args.threads)
             
        
        # Extend the assembly using the map
        if args.mode=="nanopore":
            cons_log_out= replicon.cons_log_np_out
        elif args.mode=="illumina":
            cons_log_out= replicon.cons_log_ill_out
        
        stich_telo(ref = replicon.contig_fasta,
                   left_map = replicon.l_map_out,
                   right_map = replicon.r_map_out,
                   outfile = replicon.stitch_out,
                   logout=cons_log_out,
                   tmp_left=replicon.stitch_left_fasta,
                   tmp_right=replicon.stitch_right_fasta)

    # 3: Trim consensus using a map of terminal reads onto extended
    # assembly
    # ----------------------------------------------------------------- 
    logging.info("Trimming consensus based on read support")
    
    for replicon in replicon_list:
        logging.info(f"\tContig {replicon.name}")

        if args.mode=="nanopore":
            qc_map(extended_assembly= replicon.stitch_out,
                   left = replicon.left_sam,
                   right= replicon.right_sam,
                   output_handle=replicon.trim_map,
                   t=args.threads)
            
            trim_by_map(untrimmed_assembly=replicon.stitch_out,
                        sorted_bam_file=replicon.trim_map,
                        output_handle=replicon.trim_out,
                        cons_log=cons_log_out,
                        ratio_thres=0.7,
                        qual_thres=10)
              
        elif args.mode=="illumina":
            qc_map_illumina(extended_assembly=replicon.stitch_out,
                             left_sam=replicon.left_sam,
                             right_sam=replicon.right_sam,
                             fastq_in1=args.read1,
                             fastq_in2=args.read2,
                             output_handle=replicon.trim_map,
                             t=args.threads)
            trim_by_map_illumina(untrimmed_assembly=replicon.stitch_out,
                                 sorted_bam_file = replicon.trim_map,
                                 output_handle=replicon.trim_out,
                                 cons_log = cons_log_out,
                                 cov_thres=1,
                                 ratio_thres=0.7,
                                 qual_thres=30)
    # 4: Generate QC files
    # -----------------------------------------------------------------
    logging.info("Generating QC map and finalizing result-log")

    for replicon in replicon_list:
        logging.info(f"\tContig {replicon.name}")

        if args.mode=="nanopore":
            qc_map(extended_assembly=replicon.trim_out,
                            left=replicon.left_sam,
                            right=replicon.right_sam,
                            output_handle=replicon.qc_out,
                            t=args.threads)
        if args.mode=="illumina":
            qc_map_illumina(extended_assembly=replicon.trim_out,
                             left_sam=replicon.left_sam,
                             right_sam=replicon.right_sam,
                             fastq_in1=args.read1,
                             fastq_in2=args.read2,
                             output_handle=replicon.qc_out,
                             t=args.threads)
            
        finalize_log(log = cons_log_out,
                     right_fasta = replicon.stitch_right_fasta,
                     left_fasta = replicon.stitch_left_fasta)
    
    # 5: Clean-up
    # -----------------------------------------------------------------
    logging.info("Removing temporary files")

    finished_fasta=ref_name+"_telomore.fasta"
    build_extended_fasta(org_fasta=args.reference,
                         linear_elements=linear_elements,
                         replicon_list = replicon_list,
                         output_handle=finished_fasta)
    
    shutil.move(src=finished_fasta,
                dst = os.path.join(telo_folder,finished_fasta))
    
    
    for replicon in replicon_list:
        replicon.mv_files(telo_folder,args.mode)
    
    if args.keep==False:
        
        # rmv tmp files
        for replicon in replicon_list:
            replicon.cleanup_tmp_files()     
        
        # rmv lastdb
        last_db_ext=[".bck",".des",".par",".prj",".sds",".ssp",".suf",".tis"]
        
        if args.mode=="nanopore":
            for ext in last_db_ext:
                db_file = f"{db_out}{ext}"
                os.remove(db_file)
        
        # remove map
        os.remove(map_out) # map
        os.remove(f"{map_out}.bai") # index

    logging.info(f"Output files moved to {telo_folder}")

# Run script
if __name__ == '__main__':
    
    args = get_args() # Get arguments

    setup_logging(log_file="telomore.log",quiet=args.quiet) # setup logging

    try:
        main(args)
    
    except Exception as e:
        logging.error("An error occurred during the workflow:")
        logging.error(traceback.format_exc()) 
