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

   # check if no conesnsus extens beyond the reference
   if len(left_seqs)==0:
      left_cons=SeqRecord(Seq("")) # if it is empty make an empty seqrecord to avoid errors in joining later
      print("Left cons does not extend genome")
   else:
      left_cons = SeqRecord(Seq(left_seqs[0]),id="left_cons")
      print("left cons is", len(left_cons))
   if len(right_seqs)==0:
      right_cons=SeqRecord(Seq("")) # if it is empty make an empty seqrecord to avoid errors in joining later
      print("right cons does not extend genome")
   else:
      right_cons = SeqRecord(Seq(right_seqs[0]),id="right_cons")
   
      print("right cons is", len(right_cons))
   new_genome = left_cons+genome+right_cons
   new_genome.id="Reference_with_consensus_attached"
   new_genome.description="" 
   SeqIO.write(new_genome,outfile,"fasta")
   SeqIO.write(right_cons,"tmp.right.fasta","fasta")
   SeqIO.write(left_cons,"tmp.left.fasta","fasta")
   
   # Create log of consensus length
   log = open(logout,"w")
   log.write("##############################################################################")
   log.write("\nINTIAL CONSENSUS")
   log.write("\n##############################################################################")
   log_content ="\nleft_cons:{}\tright_consensus:{}".format(len(left_cons),len(right_cons))
   log.write(log_content)
   log.write("\n>left_cons\n")
   log.write(str(left_cons.seq))
   log.write("\n>right_cons\n")
   log.write(str(right_cons.seq))
   log.close()

   return(len(left_cons),len(right_cons))

def get_support_info(bam_file, genome, position, qual_threshold=1):
   """Returns the coverage and bases matching the reference at a given position considering only base above the quality threshold"""
   fasta_file = SeqIO.read(genome,"fasta")
   bam_in = pysam.AlignmentFile(bam_file, "rb", )

   #Set read_callback="no filter" to include secondary-mappings
   # Set quality threshold=1 to include all reads
   coverage_count=bam_in.count_coverage("Reference_with_consensus_attached",start=position,stop=position+1,read_callback="nofilter",quality_threshold=1)
   A_num=coverage_count[0][0]
   C_num=coverage_count[1][0]
   T_num=coverage_count[2][0]
   G_num=coverage_count[3][0]
   cov=A_num+C_num+G_num+T_num
   
   if fasta_file.seq[position]=="N":
      matching_bases=0
   elif fasta_file.seq[position]=="A":
      matching_bases=A_num
   elif fasta_file.seq[position]=="C":
      matching_bases=C_num
   elif fasta_file.seq[position]=="G":
      matching_bases=G_num
   elif fasta_file.seq[position]=="T":
      matching_bases=T_num

   return(cov,matching_bases)

def trim_by_map(genome, sorted_bam_file, output_handle,cons_log, cov_thres=5,ratio_thres=0.7, qual_thres=0):
   """Trims the ends of a genome using an alignment of reads at the ends."""
   # load genome
   fasta = SeqIO.read(genome,"fasta")
   fasta_end = len(fasta.seq)-1 # subtract one to make it 0-indexed
   txt = open(cons_log,"r")
   txt_lines = txt.readlines()[3]
   txt.close()
   left_len = int(txt_lines.split("\t")[0].split(":")[1])
   right_len = int(txt_lines.split("\t")[1].split(":")[1])

   index_start = None
   index_end = None
   
   #trim start/left-side
   for pos in range(0,0+left_len):
      try:
         cov, match = get_support_info(sorted_bam_file,genome,pos,qual_thres)
         print("Pos",pos,"Cov",cov,"Matching",match)
         if cov>cov_thres and (match/cov)>ratio_thres:
            index_start=pos
            #print(index_start)
            break
      except TypeError: # if no reads are mapped
         continue
   
   #trim end/right
   for pos in range(fasta_end,fasta_end-right_len,-1):
      #print("POSITIONS IS ",pos)
      try:
         cov, match = get_support_info(sorted_bam_file,genome,pos,qual_thres)
         print("Pos",pos,"Cov",cov,"Matching",match)
         if cov>cov_thres and (match/cov)>ratio_thres:
            index_end=pos
            #print("INDEX END IS:",index_end)
            break
      except TypeError:
         continue
   
   # check if coverage is too low for either consensus
   # Unclear on why, but adding one on the right side is nessesary to not trim an additional base
   # Even if the consensus is rejected.
   if index_start==None and index_end==None:
      trimmed_fasta = fasta[(0+left_len):(fasta_end-right_len)+1]
      log_message = "\nLeft consensus rejected\nRight consensus rejected\n"
      trimmed_fasta.id=output_handle.split(".")[0]+"_with_no_consensus"
      trimmed_fasta.description=""
   elif index_start==None: #index without left consensus, but + right side
      log_message = "\nLeft consensus rejected\nRight consensus trimmed with {}\n".format((fasta_end-index_end))
      trimmed_fasta = fasta[(0+left_len):index_end+1]
      trimmed_fasta.id=output_handle.split(".")[0]+"_with_trimmed_consensus_attached"
      trimmed_fasta.description=""
   elif index_end==None: # index from consensus until before consensus on right side
      log_message = "\nLeft consensus trimmed with {}\nRight rejected\n".format(index_start)
      trimmed_fasta = fasta[index_start:(fasta_end-right_len)+1]
      trimmed_fasta.id=output_handle.split(".")[0]+"_with_trimmed_consensus_attached"
      trimmed_fasta.description=""
   else:
      log_message = "\nLeft consensus trimmed with {}\nRight consensus trimmed with {}\n".format(index_start,(fasta_end-index_end))
      trimmed_fasta = fasta[index_start:index_end+1]
      trimmed_fasta.id=output_handle.split(".")[0]+"_with_trimmed_consensus_attached"
      trimmed_fasta.description=""
   
   log = open(cons_log,"a")
   log.write("\n##############################################################################")
   log.write("\nCONSENSUS TRIMMING")
   log.write("\n##############################################################################")
   log.write("\nRule: Trimmed until coverage>=5 and 70% read matching position:")
   log.write(log_message)
   log.close()
   SeqIO.write(trimmed_fasta,output_handle,"fasta")

if __name__ == '__main__':

   test_map="NBC_00701.nontrimmed.map.sam.sort.bam"
   test_fasta = "NBC_00701_genome.stitch.fasta"
   trim_by_map(test_fasta,test_map,"test_out.fa","NBC_00701_consensus.log.txt")
