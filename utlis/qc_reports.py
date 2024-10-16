"""
Functions for generating useful QC metrics from the telomore script.
"""

from .cmd_tools import map_and_sort, map_and_sort_illumina
from .fasta_tools import merge_fasta, dereplicate_fastq
from .map_tools import sam_to_fastq
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import csv
import pysam

def qc_map(polished_genome,left,right,output_handle,t=1):
    '''Collect terminal reads previously identified and maps them against the extended assembly'''
    with open("all_terminal_reads.fastq", "w") as fastqfile:
        sam_to_fastq(left, fastqfile)
        sam_to_fastq(right, fastqfile)
    dereplicate_fastq("all_terminal_reads.fastq","all_terminal_reads.fastq")
    map_and_sort(polished_genome,"all_terminal_reads.fastq",output_handle,t)

def qc_map_illumina(polished_genome,left,right,output_handle,t=1):
    '''Collect terminal reads previously identified and maps them against extended assembly'''
    # Updated to sam_to_fastq instead of samtools fastq to keep
    # # non-primary reads
    with open("all_terminal_reads.fastq", "w") as fastqfile:
        sam_to_fastq(left, fastqfile)
        sam_to_fastq(right, fastqfile)
    dereplicate_fastq("all_terminal_reads.fastq","all_terminal_reads.fastq")
    map_and_sort_illumina(polished_genome,"all_terminal_reads.fastq",output_handle,t)

def cons_genome_map(left_cons,right_cons,polished_genome,output_handle,t=1):
    '''Collect consensus and maps them against the polished genome'''
    merge_fasta(left_cons,right_cons,"all_cons.fasta")
    map_and_sort(polished_genome,"all_cons.fasta",output_handle,t)

def cons_cons_map(left_cons,right_cons,output_handle,t=1):
    '''Collect consensus and maps them against each other'''
    map_and_sort(left_cons,right_cons,output_handle,t)
    
def cons_length(cons_file,output_handle,offset=100):
    '''Prints the length of fasta files'''
    cons_file = SeqIO.parse(cons_file, "fasta")
    header = ["seq_id","end_cons","full_cons"]
    tsv_log=[]
    tsv_log.append(header)

    for record in cons_file:
        seq_id = record.id
        seq_len = len(record)
        gen_len = int(seq_len)-offset
        tsv_log.append([seq_id,gen_len,seq_len])
    
    with open(output_handle, 'w', newline='') as tsv_file:
        writer = csv.writer(tsv_file,delimiter='\t')
        writer.writerows(tsv_log)

def map_to_depth(bam_file, output_handle):
    pysam.depth("-aa",bam_file,"-o",output_handle)

def finalize_log(log,right_fasta,left_fasta):
    """Function to finalize the log by prepending the final information about consensus added."""
    file = open(log)
    log_cont = file.readlines()
    file.close()
    # Org lengths of consensus added
    length_lines = log_cont[3]
    left_len = int(length_lines.split("\t")[0].split(":")[1])
    right_len = int(length_lines.split("\t")[1].split(":")[1])

    # get the number of bases trimmed off
    trim_left = log_cont[-2].split(" ")[-1]
    trim_right = log_cont[-1].split(" ")[-1]
    left_seq = SeqIO.read(left_fasta,"fasta")
    right_seq = SeqIO.read(right_fasta,"fasta")

    if trim_left.rstrip() =="rejected":
        new_left = "rejected"
        left_seq=SeqRecord(Seq(""))
    else:
        new_left = left_len - int(trim_left)
        left_seq = left_seq[int(trim_left):]
    if trim_right.rstrip() == "rejected":
        new_right="rejected"
        right_seq=SeqRecord(Seq(""))
    else:
        new_right=right_len-int(trim_right)
        right_seq = right_seq[0:new_right]

    final_lengths ="left_cons:{}\tright_consensus:{}".format(new_left,new_right)

    # write to log file
    file = open(log,"w")
    file.write("==============================================================================")
    file.write("\nFINAL GENOME EXTENSION")
    file.write("\n==============================================================================\n")
    file.write(final_lengths)
    file.write("\n>left_cons\n")
    file.write(str(left_seq.seq))
    file.write("\n>right_cons\n")
    file.write(str(right_seq.seq))
    file.write("\n")
    for line in log_cont:
        file.write(line)
    file.write("==============================================================================\n")
    file.close()
    
    #with open(filename, 'r+') as f:
    #    content = f.read()
    #    f.seek(0, 0)
    #    f.write(line.rstrip('\r\n') + '\n' + content)


if __name__ == '__main__':
   print("testing qc_report functions")
   #cons_length("all_cons.fasta","length.tsv")
   map_to_depth("np.chr.NBC_00753.fna.alignment.sam.sort.bam","teste_dreng.txt")