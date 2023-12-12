#To Do:
# Find out how to update SAM-header
# Find out whether all reads or all terminal reads is needed to train lastdb
# Maybe use the

from sam_tools import revcomp
"""
Functions for running CLI-tools, namely: minimap2, lastdb, lamassemble and medaka
"""

import subprocess

def map_and_sort(ref, fastq, t, output_handle):
   '''Maps long-reads and returns a sorted and indexed .bam-file'''
   load_command="module load minimap2/2.25"
   map_command =" ".join(["minimap2","-a","-t",str(t),str(ref),str(fastq),"-o",str(output_handle),"--secondary=yes",
                         "--secondary-seq"])
   subprocess.run(" && ".join([load_command,map_command]),shell=True, executable="/bin/bash")
   out2 = output_handle +".sort.bam"
   subprocess.run(["samtools","sort","-@",str(t),output_handle,"-o",out2])
   subprocess.run(["samtools","index",out2])

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
   
   db = db_name +".par"
   seq = open(str(output)+".fasta","w")
   subprocess.run(["lamassemble","--name="+output,db,reads],stdout=seq,)
   seq.close()

   aln = open(str(output)+".aln","w") # save alignment of reads
   subprocess.run(["lamassemble","-a",db,reads],stdout=aln)
   seq.close()


# here medaka should run
def polish(assembly, reads, output, model="r941_min_hac_g507",threads=1):
   """Polishes a fasta-file given reads and reference"""
   print("Polishing with Medaka Using", str(model))
   load_command="module load medaka"
   polish_command= " ".join(["medaka_consensus","-d",assembly,"-i",reads,"-o",output,"-m",model,"-t",str(threads)])
   subprocess.run(" && ".join([load_command,polish_command]),shell=True,executable="/bin/bash")

# if __name__ == '__main__':
#    print("Testing module functions")