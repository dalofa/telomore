"""
Functions for running CLI-tools to map reads and geneate consensus.
"""
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import subprocess
import logging
import traceback

def map_and_sort(reference: str, fastq: str, output:str,threads:int=1) ->None:
    """Maps long-reads against a reference using minimap2 through a bash
    script and returns a sorted and index bam-file"""

    # input check
    assert type(threads)==int, "threads must be an integer"
    assert os.path.isfile(reference), "the reference file specified does not exist"
    assert os.path.isfile(fastq), "the fastx-file specified does not exist"

    try:
      # run bash script
      basedir = os.path.dirname(os.path.abspath(__file__))
      # Construct the command to run the bash script using an absolute path
      cmd = " ".join(["bash", 
                     os.path.join(basedir, "..", "bash_scripts", "minimap2_cmd.sh"),
                     reference,
                     fastq,
                     str(threads),
                     output])
      minimap2_run = subprocess.run(cmd,
                                    shell=True,
                                    capture_output=True,
                                    text=True,
                                    check=True)
      # Saves the minimap2-log to the run-log
      #if minimap2_run.stderr:
         #logging.warning(f"Minimap2_log:\n{minimap2_run.stderr}")
        
    except subprocess.CalledProcessError as e:
      # If the bash script fails, capture the error and log the traceback
      logging.error(f"map_and_sort failed with error: {e}")
      logging.error(f"Script stderr: {e.stderr}")
      logging.error(traceback.format_exc())

def map_and_sort_illumina(reference, read1,read2, output,threads=1):
    """Maps illumina against a reference using bowtie2 through a bash
    script and returns a sorted and index bam-file"""
   
    # input check
    assert type(threads)==int, "threads must be an integer"
    assert os.path.isfile(reference), "the reference file specified does not exist"
    assert os.path.isfile(read1), "the fastx-file specified does not exist"
    assert os.path.isfile(read2), "the fastx-file specified does not exist"
    
    # run bash script
    basedir = os.path.dirname(os.path.abspath(__file__)) # nessesary to find the location of the bash script
    cmd = " ".join(["bash",
                    os.path.join(basedir,"..","bash_scripts","bowtie2_cmd.sh")
                    ,reference ,
                    read1,
                    read2,
                    str(threads),
                    output])
    try:
      bowtie2_run = subprocess.run(cmd,
                                   shell=True,
                                   capture_output=True,
                                   text=True,
                                   check=True)
      # Saves the bowtie2-log to the run-log
      #if bowtie2_run.stderr:
         #logging.warning(f"Minimap2_log:\n{bowtie2_run.stderr}")
      #else:
         #logging.info(bowtie2_run.stderr)

    except subprocess.CalledProcessError as e:
      # If the bash script fails, capture the error and log the traceback
      logging.error(f"map_and_sort failed with error: {e}")
      logging.error(f"Script stderr: {e.stderr}")
      logging.error(traceback.format_exc())

def map_and_sort_illumina_cons(reference, consensus_fasta, output,threads=1):
    """Maps consensus against a reference using bowtie2 through a bash
    script and returns a sorted and index bam-file"""

    # input check
    assert type(threads)==int, "threads must be an integer"
    assert os.path.isfile(reference), "the reference file specified does not exist"
    assert os.path.isfile(consensus_fasta), "the fastx-file specified does not exist"

    try:
      # run bash script
      basedir = os.path.dirname(os.path.abspath(__file__))
      # Construct the command to run the bash script using an absolute path
      cmd = " ".join(["bash", 
                     os.path.join(basedir, "..", "bash_scripts", "bowtie2_cons_cmd.sh"),
                     reference,
                     consensus_fasta,
                     str(threads),
                     output])
      bowtie2_cons_run = subprocess.run(cmd,
                                    shell=True,
                                    capture_output=True,
                                    text=True,
                                    check=True)

    except subprocess.CalledProcessError as e:
      # If the bash script fails, capture the error and log the traceback
      logging.error(f"map_and_sort_illumina_cons failed with error: {e}")
      logging.error(f"Script stderr: {e.stderr}")
      logging.error(traceback.format_exc())
      

def train_lastDB(fasta_name,reads,db_name, t=1):
   '''Trains and lastDB database using a reference and long-reads'''
   # index fasta file
   try:
      lastdb_run = subprocess.run(["lastdb",
                                   "-P"+str(t),
                                   "-uRY4",
                                   db_name,
                                   fasta_name],
                                   check=True)
   except subprocess.CalledProcessError as e:
      # If the bash script fails, capture the error and log the traceback
      logging.error(f"train_lastDB failed with error: {e}")
      logging.error(f"Script stderr: {e.stderr}")
      logging.error(traceback.format_exc())
   
   # train last-db
   try:
      file = open(db_name +".par","w") # needed to write to file
      subprocess.run(["last-train",
                      "-P"+str(t),
                      "-Qkeep",
                      db_name,
                      reads],
                      stdout=file,
                      check=True)
      file.close()
   except subprocess.CalledProcessError as e:
      # If the bash script fails, capture the error and log the traceback
      logging.error(f"train_lastDB failed with error: {e}")
      logging.error(f"Script stderr: {e.stderr}")
      logging.error(traceback.format_exc())

def generate_consensus_lamassemble(db_name,reads,output,flip="no"):
   """Generates a consensus fasta-file given a LAST-DB and a series of reads. 
   Is suitable for dissimilar reads such as nanopore"""
   
   # Check if only a single read is present before generating consensus
   # If that is the case write just that read to consensus file
   # This prevents issues down the line when stitching chromosome and 
   sequence_count = 0
   for record in SeqIO.parse(reads, "fastq"):
      sequence_count += 1

   if sequence_count==0:
      logging.info(f"There are no reads to construct a consensus from. Emtpy consensus returned to {output}")
      with open(f"{output}", "w") as seq:
         empty_record = SeqRecord(Seq(""), id="empty_consensus")
         SeqIO.write(empty_record, seq, "fasta")
   if sequence_count == 1:
      single_record = record
      logging.info(f"There are only a single read to construct a consensus from. Returning read as consensus to {output}")
      with open(f"{output}", "w") as seq:
         SeqIO.write(single_record, seq, "fasta")
   elif sequence_count>1:
      db = db_name +".par"
      seq = open(str(output),"w")
      try:
         lamassemble_run = subprocess.run(["lamassemble",
                                           "--name="+output,
                                           db,
                                           reads],
                                           stdout=seq,
                                           check=True)
         seq.close()

         aln = open(str(output)+".aln","w") # save alignment of reads
         subprocess.run(["lamassemble","-a",db,reads],stdout=aln)
         aln.close()

      except subprocess.CalledProcessError as e:
      # If the bash script fails, capture the error and log the traceback
         logging.error(f"generate_consensus_lamassemble failed with error: {e}")
         logging.error(f"Script stderr: {e.stderr}")
         logging.error(traceback.format_exc())

def generate_consensus_mafft(reads, output):
   """Generates a consensus fasta-file given a series of reads in fasta-format. 
   Relies on very similar reads."""

   #Check if only a single read is mapped and use that as the consensus if there are no others.
   sequence_count = 0
   for record in SeqIO.parse(reads, "fastq"):
      sequence_count += 1

   if sequence_count==0:
      logging.info(f"There are no reads to construct a consensus from. Emtpy consensus returned to {output}")
      with open(f"{output}", "w") as seq:
         empty_record = SeqRecord(Seq(""), id="empty_consensus")
         SeqIO.write(empty_record, seq, "fasta")
   
   if sequence_count == 1:
      single_record = record
      logging.info(f"There are only a single read to construct a consensus from. Returning read as consensus to {output}")
      with open(f"{output}", "w") as seq:
         SeqIO.write(single_record, seq, "fasta")
   
   elif sequence_count>1:
      # Convert fastq to fasta
      fasta_reads = reads +".fasta"
      with open(reads, "r") as input_handle, open(fasta_reads, "w") as output_handle:
         SeqIO.convert(input_handle, "fastq", output_handle, "fasta")
      # Run Mafft
      mafft_output=reads+".aln"
      try:
         mafft_file=open(mafft_output,"w")
         mafft_run = subprocess.run(["mafft","--quiet",
                                  fasta_reads],
                                  stdout=mafft_file,
                                  check=True)
         mafft_file.close()
      
      except subprocess.CalledProcessError as e:
      # If the bash script fails, capture the error and log the traceback
         logging.error(f"generate_consensus_mafft failed with error: {e}")
         logging.error(f"Script stderr: {e.stderr}")
         logging.error(traceback.format_exc())


      try:
            # Generate consensus using Emboss cons
         # Plurality 1 ensures that only one reads needs to cover a position to generate consensus
         # This is dangerous with dissimilar seqeunces
         subprocess.run(["cons",
                         "-name="+output
                         ,"-plurality",str(1),
                         "-sequence="+mafft_output,
                         "-outseq="+output],
                         capture_output=True,
                         text=True,
                         check=True)
      except subprocess.CalledProcessError as e:
      # If the bash script fails, capture the error and log the traceback
         logging.error(f"generate_consensus_mafft failed with error: {e}")
         logging.error(f"Script stderr: {e.stderr}")
         logging.error(traceback.format_exc())

if __name__ == '__main__':
   print("Testing module functions")

   map_and_sort_illumina_cons(reference="NBC_00357_telomore_extended.chrom.left.fa",
                              consensus_fasta="left_cons.fasta", 
                              threads=16,
                              output="blahhh")
