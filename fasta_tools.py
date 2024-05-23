from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def get_chromosome(fasta, output_handle):
    '''Returns the longest contig in a fasta-file'''

    # test if there are a single entry in the fasta file
    try: # there is a single entry
        chromosome = SeqIO.read(fasta,format="fasta")
        SeqIO.write(chromosome,output_handle, format="fasta")
        message = "A single contig: {} was found and will be used for mapping".format(">"+chromosome.id)
        print(message)

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
        print(message)


def attach_seq(left,right,chromsome,output_name,offset=0):
    '''Attaches a left and right fasta sequence to a cromosome'''

    l = SeqIO.read(left,"fasta")
    r = SeqIO.read(right,"fasta")
    chrom = SeqIO.read(chromsome,"fasta")
    
    if offset==0: # if offset is 0 offset:-offset fucks it up
        genome = chrom
    elif offset >=len(chrom.seq)/2:
        print("Error: Offset is larger than  1/2 genome length.")
        return
    else:
        genome = chrom[offset:-offset]

    att_genome = l + genome + r
    att_genome.id=output_name.split(".")[0]
    SeqIO.write(att_genome,output_name, "fasta")

# A function to merge fasta files
def merge_fasta(input_file1, input_file2, output_file):
    '''Merges two fasta files into one'''
    # Read sequences from input_file1 and input_file2
    sequences1 = list(SeqIO.parse(input_file1, "fasta"))
    sequences2 = list(SeqIO.parse(input_file2, "fasta"))

    # Merge sequences
    merged_sequences = sequences1 + sequences2

    # Write merged sequences to output_file
    with open(output_file, "w") as output_handle:
        SeqIO.write(merged_sequences, output_handle, "fasta")

### ADD OPTION TO TRIM FROM START OR END

def trim_to_cons(input_seq,num_base,output_handle):
    "Trims fasta file to the first x_num of bases in it"
    
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


                #print(seq[0:num_base+1])
            else:
                print("Error: Index out of range")
        
        if len(all_rec)>0:
            SeqIO.write(all_rec,output_handle,"fasta")



def strip_fasta(input_file, output_file, x, remove_from='start'):
    """
    Modify sequences in a FASTA file by removing the first or last x bases.
    """

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



if __name__ == '__main__':
    import os
    #print(os.getcwd())
    #attach_seq("left.fasta","right.fasta","middle.fasta","new3.fasta",offset=40)
    trim_to_cons("3_cons.fasta",120,"ot.boi")
    