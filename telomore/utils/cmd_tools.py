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
from pkg_resources import resource_filename


def map_and_sort(reference: str, fastq: str, output:str,threads:int=1) -> None:
    """Maps long-reads against a reference using minimap2 through a bash
    script and returns a sorted and index bam-file"""
    # input check
    assert type(threads)==int, "threads must be an integer"
    assert os.path.isfile(reference), "the reference file specified does not exist"
    assert os.path.isfile(fastq), "the fastx-file specified does not exist"

    try:
      # run bash script
      # Use resource_filename to make compatible with python packaging
      script_path = resource_filename('telomore', 'bash_scripts/minimap2_cmd.sh')
      cmd = " ".join(["bash", 
                     script_path,
                     reference,
                     fastq,
                     str(threads),
                     output])
      minimap2_run = subprocess.run(cmd,
                                    shell=True,
                                    capture_output=True,
                                    text=True,
                                    check=True)
    except subprocess.CalledProcessError as e:
      # If the bash script fails, capture the error and log the traceback
      logging.error(f"map_and_sort failed with error: {e}")
      logging.error(f"Script stderr: {e.stderr}")
      logging.error(traceback.format_exc())
      raise

def map_and_sort_illumina(reference: str, read1: str, read2: str, output: str,threads=1) -> None:
    """Maps illumina against a reference using bowtie2 through a bash
    script and returns a sorted and index bam-file"""
   
    # input check
    assert type(threads)==int, "threads must be an integer"
    assert os.path.isfile(reference), "the reference file specified does not exist"
    assert os.path.isfile(read1), "the fastx-file specified does not exist"
    assert os.path.isfile(read2), "the fastx-file specified does not exist"
    
    # run bash script
    script_path = resource_filename('telomore', 'bash_scripts/bowtie2_cmd.sh')
    cmd = " ".join(["bash",
                    script_path
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
    except subprocess.CalledProcessError as e:
      # If the bash script fails, capture the error and log the traceback
      logging.error(f"map_and_sort failed with error: {e}")
      logging.error(f"Script stderr: {e.stderr}")
      logging.error(traceback.format_exc())
      raise

def map_and_sort_illumina_cons(reference:str, consensus_fasta:str, output:str,threads=1) -> None:
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
      raise

def train_lastDB(fasta_name:str,reads:str,db_name:str, t=1) -> None:
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
      raise
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
      raise

def generate_consensus_lamassemble(db_name:str,reads:str,output:str) -> None:
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
         raise

def generate_consensus_mafft(reads:str, output:str) -> None:
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
      mafft_output=output+".aln"
      try:
         mafft_file=open(mafft_output,"w")
         mafft_run = subprocess.run(["mafft","--quiet",
                                  fasta_reads],
                                  stdout=mafft_file,
                                  check=True)
         mafft_file.close()
         os.remove(fasta_reads)
      
      except subprocess.CalledProcessError as e:
      # If the bash script fails, capture the error and log the traceback
         logging.error(f"generate_consensus_mafft failed with error: {e}")
         logging.error(f"Script stderr: {e.stderr}")
         logging.error(traceback.format_exc())
         raise

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
         raise
