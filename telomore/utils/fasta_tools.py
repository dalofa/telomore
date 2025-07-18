"""
Utilities for handling fasta files.
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging

def get_linear_elements(fasta_file:str) -> list:
    """Returns a list of contigs that have linear in their header description"""

    linear_list = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        if "linear" in record.description:
            linear_list.append(record.id)
    return linear_list
    
def extract_contig(fasta_in:str,contig_name:str,fasta_out:str):
    """Writes a fasta file of a specified contig"""
    for record in SeqIO.parse(fasta_in, "fasta"):
        if record.id==contig_name:
            contig = record
            with open(fasta_out,"w") as fq_file:
                SeqIO.write(sequences= contig,
                            handle= fasta_out,
                            format="fasta")

def get_fasta_length(fasta_file:str,contig_name:str) -> int:
    """Returns the length of a fasta record specified by the user"""
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id==contig_name:
            length = len(record.seq)
            return(length)

def dereplicate_fastq(fastq_in:str, fastq_out:str) -> None:
    """Dereplicates a fastq file by read-name and sequence."""
    seen_reads = set()  # To store unique read identifiers and sequences
    unique_reads = []

    with open(fastq_in, "r") as infile:
        for record in SeqIO.parse(infile, "fastq"):
            # Combine read ID and sequence to ensure uniqueness by both criteria
            read_key = (record.id, str(record.seq))
            
            if read_key not in seen_reads:
                seen_reads.add(read_key)
                unique_reads.append(record)

    with open(fastq_out, "w") as outfile:
        SeqIO.write(unique_reads, outfile, "fastq")

def cat_and_derep_fastq(fastq_in1:str,fastq_in2:str,fastq_out:str) -> None:
    """Concatanates two fastq-files and dereplicates it"""
    with open(fastq_out, "w") as outfile:
    
        # concat 
        with open(fastq_in1, "r") as infile1:
            for record in SeqIO.parse(infile1, "fastq"):
                SeqIO.write(record, outfile, "fastq")
    
        with open(fastq_in2, "r") as infile2:
            for record in SeqIO.parse(infile2, "fastq"):
                SeqIO.write(record, outfile, "fastq")

    dereplicate_fastq(fastq_in = fastq_out,
                      fastq_out=fastq_out)


def get_chromosome(fasta:str, output_handle:str) -> None:
    """Writes the longest contig to a fasta-file"""

    # test if there are a single entry in the fasta file
    try: # there is a single entry
        chromosome = SeqIO.read(fasta,format="fasta")
        SeqIO.write(chromosome,output_handle, format="fasta")
        message = "A single contig: {} was found and will be used for mapping".format(">"+chromosome.id)
        logging.info(message)

    except ValueError: # there are more than one entry
        contigs = SeqIO.parse(fasta,format="fasta")
        max_len = 0

        # identify longest entry and assume it is the chromosome
        for record in contigs:
            seq_len = len(record.seq)
            
            if seq_len > max_len:
                chromosome = record
                max_len = seq_len
        
        SeqIO.write(chromosome,output_handle, format="fasta")
        message = "The longest contig: {} has been saved as {} and will be used for mapping.".format(">"+chromosome.id,output_handle)
        logging.info(message)


def attach_seq(left:str,right:str,chromsome:str,output_name:str,offset: int=0) -> None:
    """Attaches a left and right fasta sequence to an assembly sequence and writes it to a new file."""

    l = SeqIO.read(left,"fasta")
    r = SeqIO.read(right,"fasta")
    chrom = SeqIO.read(chromsome,"fasta")
    
    if offset==0: # if offset is 0 offset:-offset fucks it up
        genome = chrom
    elif offset >=len(chrom.seq)/2:
        logging.error("Error: Offset is larger than  1/2 genome length.")
        return
    else:
        genome = chrom[offset:-offset]

    att_genome = l + genome + r
    att_genome.id=output_name.split(".")[0]
    SeqIO.write(att_genome,output_name, "fasta")

# A function to merge fasta files
def merge_fasta(input_file1:str, input_file2:str, output_file:str) -> None:
    """Merges two fasta files into one"""
    # Read sequences from input_file1 and input_file2
    sequences1 = list(SeqIO.parse(input_file1, "fasta"))
    sequences2 = list(SeqIO.parse(input_file2, "fasta"))

    # Merge sequences
    merged_sequences = sequences1 + sequences2

    # Write merged sequences to output_file
    with open(output_file, "w") as output_handle:
        SeqIO.write(merged_sequences, output_handle, "fasta")

def trim_to_cons(input_seq:str,num_base:int,output_handle:str) -> None:
    """Trims fasta file to the first x_num of bases in it and writes it to a new file."""
    
    # load file
    with open(input_seq) as fasta_file:
        
        all_rec = []

        for record in SeqIO.parse(fasta_file,"fasta"):
            new_id = "trimmed_"+record.id

            r_seq = record.seq
            length = len(r_seq)

            if num_base <= length:
                to_write = SeqRecord(seq=r_seq[0:num_base+1],id=new_id,description="")
                all_rec.append(to_write)

            else:
                logging.error("Error: Index out of range")
        
        if len(all_rec)>0:
            SeqIO.write(all_rec,output_handle,"fasta")



def strip_fasta(input_file:str, output_file:str, x:int, remove_from:str='start') -> None:
    """Strip fasta-file of the first of last x-bases."""
    
    assert type(x)==int
    
    records = []
    
    for record in SeqIO.parse(input_file, "fasta"):
        if remove_from == 'start':
            modified_seq = record.seq[x:]
        elif remove_from == 'end':
            modified_seq = record.seq[:-x]
        else:
            raise ValueError("remove_from must be either 'start' or 'end'")
        
        record.seq = modified_seq
        records.append(record)
    
    SeqIO.write(records, output_file, "fasta")


def build_extended_fasta(org_fasta:str, linear_elements: list, replicon_list:list,output_handle:str) -> None:
    """Rejoins extended contigs in the order of the original fasta and writes it to a new file."""
    
    seq_rec_list = [] # list of seqrecord to write to newfile

    for record in SeqIO.parse(org_fasta, "fasta"):
        if record.id in linear_elements:
            for replicon in replicon_list:
                if replicon.name==record.id:
                    path_to_telomore_rec = replicon.trim_out
            telomore_rec = SeqIO.read(path_to_telomore_rec,format="fasta")
            telomore_rec.description= "[linear]"
            seq_rec_list.append(telomore_rec)
        else:
            seq_rec_list.append(record)
    
    SeqIO.write(sequences=seq_rec_list,
                handle=output_handle,
                format="fasta")
