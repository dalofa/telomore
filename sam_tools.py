#To-Do: Find a way to update sam-headers, such that it is clear which filteres have been applied

import pysam
import re
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

###_______Functions________
def get_terminal_reads(sam_file,loutput_handle,routput_handle):

   """A function that retrieves all reads mapping at the very start or end of a reference"""
   
   input = pysam.AlignmentFile(sam_file, "r")
   
   # Fetch all reads aligned at start or end of reference
   seq_end = input.lengths[0]
   ref_name = input.get_reference_name(0)
   left_reads = input.fetch(ref_name,start=0, stop=1)
   right_reads = input.fetch(ref_name,seq_end-1,seq_end)

   # Write all fetched reads to a new file
   lterminal_reads = pysam.AlignmentFile(loutput_handle,"w",template=input)
   rterminal_reads = pysam.AlignmentFile(routput_handle,"w",template=input)

   for lread in left_reads:
       if lread.query_sequence is None:
          continue
       else:
          lterminal_reads.write(lread)

   for rread in right_reads:
       if rread.query_sequence is None:
          continue
       else:
          rterminal_reads.write(rread)

   lterminal_reads.close()
   rterminal_reads.close()
   input.close()

def get_left_soft(sam_file,left_out,offset=0,):
    """A function that retrieves reads soft-clipped at 5'-end
    It returns the reads as a sam file and only soft-clipped part as fastq-file"""
    sam_in = pysam.AlignmentFile(sam_file, "r")
    lclip = pysam.AlignmentFile(left_out+".sam","w",template=sam_in)
    lfastq = open(left_out+".fastq","w")
    
    start_clip = r"^(\d+)S"
    for read in sam_in:
        lmatch = re.match(start_clip,read.cigarstring)
        
        if lmatch:
           clip_num=int(lmatch.group(1)) # digits are retrieve via .group

           if clip_num >read.reference_start:
              lclip.write(read) # write to sam-file

              # get info for fastq-file
              name = read.query_name
              seq = read.query_sequence[0:clip_num+offset]
              sanger_qual = "".join([chr(q+33) for q in read.query_qualities[0:clip_num+offset]]) # phred qual converted to ASCII with 33 offset
              lfastq.write("@{}\n{}\n+\n{}\n".format(name, seq, sanger_qual))
    sam_in.close()
    lclip.close()
    lfastq.close()

def get_right_soft(sam_file,left_out,offset=0,):
    """A function that retrieves reads soft-clipped at 3'-end
    It returns the reads as a sam file and only soft-clipped part as fastq-file"""
    sam_in = pysam.AlignmentFile(sam_file, "r")
    rclip = pysam.AlignmentFile(left_out+".sam","w",template=sam_in)
    rfastq = open(left_out+".fastq","w")
    seq_end = sam_in.lengths[0] # get length of reference
    end_clip = r"(\d+)S$"
    for read in sam_in:
        rmatch = re.search(end_clip,read.cigarstring)
        if rmatch:
           clip_num=int(rmatch.group(1)) # digits are retrieve via .group
           
           if clip_num+read.reference_end>seq_end:
              rclip.write(read) # write to sam-file
              
              # get info for fastq-file
              name = read.query_name
              seq = read.query_sequence[-(clip_num+offset):]
              sanger_qual = "".join([chr(q+33) for q in read.query_qualities[-(clip_num+offset):]]) # phred qual converted to ASCII with 33 offset
              rfastq.write("@{}\n{}\n+\n{}\n".format(name, seq, sanger_qual))

    sam_in.close()
    rclip.close()
    rfastq.close()

def revcomp_reads(reads_in,reads_out):
   """A function that takes in fastq reads and retunrs their reverse complement"""
   
   with open(reads_in, "r") as input_handle, open(reads_out, "w") as output_handle:
    
    for record in SeqIO.parse(input_handle, "fastq"):
        # Get the reverse complement of the sequence
        rev_complement_seq = record.seq.reverse_complement()
        
        # Reverse the quality scores as well
        rev_quality_scores = record.letter_annotations["phred_quality"][::-1]
        
        # Create a new record with the reverse complement sequence and quality scores
        rev_complement_record = record
        rev_complement_record.id = "rev_"+str(record.id)


        rev_complement_record.seq = rev_complement_seq
        rev_complement_record.letter_annotations["phred_quality"] = rev_quality_scores
        
        # Write the reverse complement record to the output FASTQ file
        SeqIO.write(rev_complement_record, output_handle, "fastq")
        
def revcomp(fasta_in,fasta_out):
   with open(fasta_in, "r") as input_handle, open(fasta_out, "w") as output_handle:

      for record in SeqIO.parse(input_handle, "fasta"):
        # Get the reverse complement of the sequence
        rev_complement_seq = record.seq.reverse_complement()
        
        # Create a new record with the reverse complement sequence and quality scores
        rev_complement_record = record
        rev_complement_record.id = "rev_"+str(record.id)
        rev_complement_record.seq = rev_complement_seq
        
        # Write the reverse complement record to the output FASTQ file
        SeqIO.write(rev_complement_record, output_handle, "fasta")

def stich_telo(ref,left_map,right_map,outfile,logout="consensus.log.txt"):
   # function to stitch telomeres onto the fucking reference
   # add checks for mapping positions

   # extract left cons-to-stich
   l_sam_in = pysam.AlignmentFile(left_map, "r")
   left_seqs =[]
   start_clip = r"^(\d+)S"

   # filter away mapping at right side
   cons_at_left = [read for read in l_sam_in if read.reference_start<1000]

   for read in cons_at_left:
      lmatch = re.match(start_clip,read.cigarstring)
      if lmatch:
         clip_num=int(lmatch.group(1)) # digits are retrieve via .group
         seq = read.query_sequence[0:(clip_num-read.reference_start)] # Adjust for if more than just overhanging bases are soft-clipped
         left_seqs.append(seq)
   l_sam_in.close()

   # extract right cons-to-stich
   r_sam_in = pysam.AlignmentFile(right_map, "r")
   seq_end = r_sam_in.lengths[0] # get length of reference
   right_seqs = []
   end_clip = r"(\d+)S$" # reg. exp for ending with *S[num]

   cons_at_right = [read for read in r_sam_in if read.reference_start>1000]
   for read in cons_at_right:
      rmatch = re.search(end_clip,read.cigarstring)
      if rmatch:

         clip_num=int(rmatch.group(1)) # digits are retrieve via .group
         # Adjusting for potential difference between overhang and soft-clip
         adj = seq_end-read.reference_end
         if clip_num+read.reference_end>seq_end:
            seq = read.query_sequence[-(clip_num-adj):]
            right_seqs.append(seq)
   r_sam_in.close()

   # stich the fuckers toghether
   genome = SeqIO.read(ref,"fasta")
   left_cons = SeqRecord(Seq(left_seqs[0]),id="left_cons")
   right_cons = SeqRecord(Seq(right_seqs[0]),id="right_cons")
   print("left cons is", len(left_cons))
   print("right cons is", len(right_cons))
   new_genome = left_cons+genome+right_cons
   new_genome.id="Reference_with_consensus_attached"
   new_genome.description="" 
   SeqIO.write(new_genome,outfile,"fasta")
   
   # Create log of consensus length
   log = open(logout,"w")
   log_content ="left_cons:{}\tright_consensus:{}".format(len(left_cons),len(right_cons))
   log.write(log_content)
   log.write("\n>left_cons\n")
   log.write(left_cons.seq)
   log.write("\n>right_cons\n")
   log.write(right_cons.seq)
   log.close()
   
   return(len(left_cons),len(right_cons))

def get_support_info(bam_file, genome, position, qual_threshold=1):
   """Returns the coverage and bases matching the reference at a given position considering only base above the quality threshold"""

   bam_in = pysam.AlignmentFile(bam_file, "rb", )
   
   # pileup: Truncate=True ensures that only the specific position supplied is returned and not full length reds
   # stepper="nofilter" ensures that secondary alignments are also returned
   pileup = bam_in.pileup( start=position, stop=position+1,truncate=True,stepper="nofilter")
   fasta_file = SeqIO.read(genome,"fasta")
   # counter
   matching_bases = 0

   for pileup_column in pileup:
      pileup_column.set_min_base_quality(qual_threshold) # must be applied to not use the default base qual filter which is >13
      cov = len(pileup_column.get_query_qualities()) # .get_query_qualities/get_query_sequences returns only reads passsing the filters
      pos_base = pileup_column.get_query_sequences()
      #print(pos_base)
      # check how many bases match the reference
      for base in pos_base:
         if fasta_file.seq[position].upper() == base.upper():
            matching_bases+=1
      return(cov,matching_bases)


def trim_by_map(genome, sorted_bam_file, output_handle,cons_log, cov_thres=5,ratio_thres=0.7, qual_thres=1):
   """Trims the ends of a genome using an alignment of reads at the ends."""
   # load genome
   fasta = SeqIO.read(genome,"fasta")
   fasta_end = len(fasta.seq)-1
   txt = open(cons_log,"r")
   txt_lines = txt.readlines()[0]
   txt.close()
   left_len = int(txt_lines.split("\t")[0].split(":")[1])
   right_len = int(txt_lines.split("\t")[1].split(":")[1])

   index_start = None
   index_end = None
   
   #trim start/left-side
   for pos in range(0,0+left_len):
      try:
         cov, match = get_support_info(sorted_bam_file,genome,pos,qual_thres)
         if cov>cov_thres and (match/cov)>ratio_thres:
            index_start=pos
            break
      except TypeError: # if no reads are mapped
         continue
   
   #trim end/right
   for pos in range(fasta_end,fasta_end-right_len,-1):
      try:
         cov, match = get_support_info(sorted_bam_file,genome,pos,qual_thres)
         if cov>cov_thres and (match/cov)>ratio_thres:
            index_end=pos
            break
      except TypeError:
         continue
   
   # check if coverage is too low for either consensus
   if index_start==None and index_end==None:
      trimmed_fasta = fasta[(0+left_len):(fasta_end-right_len)]
      print("Coverage too low: both consensus rejected")
      log_message = "\nBoth consensus rejceted and genome wihtout extension returned."
   elif index_start==None:
      print("Coverage too low for left consensus")
      #index without consensus + right side
      log_message = "\nLeft consensus rejected. Right consensus trimmed with {}".format(fasta_end-index_end)
   elif index_end==None:
      print("Coverage too low for right consensus")
      log_message = "\nRight consensus rejected. Left consensus trimmed with {}".format(index_start)
      # index from consensus until before consensus on right side
   else:
      log_message = "\nLeft consensus trimmed with {}.\nRight consensus trimmed with {}".format(index_start,(fasta_end-index_end))
      trimmed_fasta = fasta[index_start:index_end]
      trimmed_fasta.id=output_handle.split(".")[0]+"_with_trimmed_consensus_attached"
      trimmed_fasta.description=""
   
   txt = open(cons_log,"a")
   txt.write(log_message)
   txt.close()
   SeqIO.write(trimmed_fasta,output_handle,"fasta")

if __name__ == '__main__':
   print("testing samfilter functions")
   ref = "Reference_with_consensus_attached"
   test_file="NBC_00287.nonpolish.qc.map.sam.sort.bam"
   genome="NBC_00287_genome.stitch.fasta"
   test_output="test.fa"
   cons_log="consensus.log.txt"
   trim_by_map(genome,test_file,test_output,cons_log="consensus.log.txt")

