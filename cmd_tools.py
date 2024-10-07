#To Do:
# Find out how to update SAM-header
# Find out whether all reads or all terminal reads is needed to train lastdb
# Maybe use the

from sam_tools import revcomp
import os
from Bio import SeqIO


"""
Functions for running CLI-tools, namely: minimap2, lastdb, lamassemble and medaka
"""

import subprocess

def map_and_sort(reference, fastq, output,threads=1):
    """Maps long-reads against a reference using minimap2 through a bash
    script and returns a sorted and index bam-file"""

    # input check
    assert type(threads)==int, "threads must be an integer"
    assert os.path.isfile(reference), "the reference file specified does not exist"
    assert os.path.isfile(fastq), "the fastx-file specified does not exist"
    # run bash script
    basedir = os.path.dirname(__file__) # nessesary to find the location of the bash script
    cmd = " ".join(["bash", os.path.join(basedir,"minimap2_cmd.sh") ,reference , fastq , str(threads), output])
    subprocess.run(cmd,shell=True)

def map_and_sort_illumina(reference, fastq, output,threads=1):
    """Maps illumina against a reference using bowtie2 through a bash
    script and returns a sorted and index bam-file"""
    # input check
    assert type(threads)==int, "threads must be an integer"
    assert os.path.isfile(reference), "the reference file specified does not exist"
    assert os.path.isfile(fastq), "the fastx-file specified does not exist"
    
    # run bash script
    basedir = os.path.dirname(__file__) # nessesary to find the location of the bash script
    cmd = " ".join(["bash", os.path.join(basedir,"bowtie2_cmd.sh") ,reference , fastq , str(threads), output])
    subprocess.run(cmd,shell=True)
    

def train_lastDB(fasta_name,reads,db_name, t=1):
   '''Trains and lastDB database using a reference and long-reads'''
   # index fasta file
   print("lastdb: Indexing genome")
   subprocess.run(["lastdb","-P"+str(t),"-uRY4",db_name,fasta_name])
   print("last-train: Determining substituion and gap rates")
   # train last-db
   file = open(db_name +".par","w") # needed to write to file
   subprocess.run(["last-train","-P"+str(t),"-Qkeep",db_name,reads], stdout=file)
   file.close()

def generate_consensus(db_name,reads,output,flip="no"):
   """Generates a consensus fasta-file given a LAST-DB and a series of reads."""
   
   # Check if only a single read is present before generating consensus
   # If that is the case write just that read to consensus file
   # This prevents issues down the line when stitching chromosome and 
   sequence_count = 0
   for record in SeqIO.parse(reads, "fastq"):
      sequence_count += 1
   if sequence_count == 1:
      single_record = record
      print("There are only a single read to construct a consensus from. Returning read as consensus.")
      with open(f"{output}.fasta", "w") as seq:
         SeqIO.write(single_record, seq, "fasta")  
   elif sequence_count>1:

      db = db_name +".par"
      seq = open(str(output)+".fasta","w")
      subprocess.run(["lamassemble","--name="+output,db,reads],stdout=seq,)
      seq.close()

      aln = open(str(output)+".aln","w") # save alignment of reads
      subprocess.run(["lamassemble","-a",db,reads],stdout=aln)
      aln.close()


# here medaka should run
def polish(assembly, reads, output, model="r941_min_hac_g507",threads=1):
   """Polishes a fasta-file given reads and reference"""
   print("Polishing with Medaka Using", str(model))
   load_command="module load medaka"
   polish_command= " ".join(["medaka_consensus","-d",assembly,"-i",reads,"-o",output,"-m",model,"-t",str(threads)])
   subprocess.run(" && ".join([load_command,polish_command]),shell=True,executable="/bin/bash")

# if __name__ == '__main__':
#    print("Testing module functions")