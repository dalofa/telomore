from cmd_tools import map_and_sort
from fasta_tools import merge_fasta
from Bio import SeqIO
import subprocess
import csv
import pysam

"""
Functions for generating useful QC metrics from the telomore script.
"""


def qc_map(polished_genome,left,right,output_handle,t=1):
   '''Collect terminal reads previously identified and maps them against the polished genome'''
   collect_fastq = open("all_terminal_reads.fastq","a")
   subprocess.run(["samtools","fastq","-@",str(t),left],stdout=collect_fastq)
   subprocess.run(["samtools","fastq","-@",str(t),right],stdout=collect_fastq)
   collect_fastq.close()
   map_and_sort(polished_genome,"all_terminal_reads.fastq",output_handle,t)

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



if __name__ == '__main__':
   print("testing qc_report functions")
   #cons_length("all_cons.fasta","length.tsv")
   map_to_depth("np.chr.NBC_00753.fna.alignment.sam.sort.bam","teste_dreng.txt")