"""
Functions for handling read mappings. Primarily for extracting and filtering terminal reads.
"""

import pysam
import re
import gzip
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def map_to_contig(bam_in,contig_name):
   """Takes a .bam file for multiple contigs and extracts the mapping for the contig"""
   return

def sam_to_readpair(sam_in,fastq_in1, fastq_in2,fastq_out1,fastq_out2):
   with pysam.AlignmentFile(sam_in) as samfile:
      reads_to_grep = set() # using a set should be faster than list

      # get all read names
      for read in samfile.fetch(until_eof=True):
         read_name = read.query_name 
         if " " in read.query_name:
            read_name = read_name.split(" ")[0]
         
         reads_to_grep.add(read_name)
      
      # get read 1
      with gzip.open(fastq_in1,"rt") as gzip_handle, open(fastq_out1,"w") as outfile:
         for record in SeqIO.parse(gzip_handle,"fastq"):
            if record.id in reads_to_grep:
               SeqIO.write(record, outfile, "fastq")

      # get read 2
      with gzip.open(fastq_in2,"rt") as gzip_handle, open(fastq_out2,"w") as outfile:
         for record in SeqIO.parse(gzip_handle,"fastq"):
            if record.id in reads_to_grep:
               SeqIO.write(record, outfile, "fastq")

def sam_to_fastq(sam_in,fastq_out):
   """Convert a sam-file to fastq-format, excluding unmapped reads."""
   with pysam.AlignmentFile(sam_in, "r") as samfile:
         for read in samfile.fetch(until_eof=True):
               if not read.is_unmapped:
                  name = read.query_name
                  seq = read.query_sequence
                  qual = read.qual
                  if qual is None:
                     qual = 'I' * len(seq)  # Assign a default high-quality score
                  
                  # Write the read in FASTQ format to the provided handle
                  fastq_out.write(f"@{name}\n{seq}\n+\n{qual}\n")

def mapped_bases(cigarstring):
   """Calculate the number of bases mapped to the reference in a given CIGAR string."""
   # Define operations that consume reference bases
   consuming_operations = "MDNX="

   # Parse the CIGAR string using regex
   # This produces a tuple in the format (121,"S")
   operations = re.findall(r'(\d+)([MIDNSHP=X])', cigarstring)
   
   # Initialize base count
   mapped_bases_count = 0

   # Loop through the parsed operations and sum bases for consuming operations
   for length, op in operations:
      if op in consuming_operations:
            mapped_bases_count += int(length)
   
   return mapped_bases_count

def cigar_maps_more_bases(cigar1, cigar2):
   """Compare two CIGAR strings and determine which one maps to more bases."""
   bases1 = mapped_bases(cigar1)
   bases2 = mapped_bases(cigar2)
   
   if bases1 > bases2:
      return True
   elif bases1 < bases2:
      return False

def get_terminal_reads(sam_file,loutput_handle,routput_handle):
   """A function that retrieves all reads mapping at the very start or end of a reference"""
   
   input = pysam.AlignmentFile(sam_file, "r")
   
   # Fetch all reads aligned at start or end of reference
   seq_end = input.lengths[0]
   ref_name = input.get_reference_name(0)
   left_reads = input.fetch(ref_name,start=0, stop=20)
   right_reads = input.fetch(ref_name,seq_end-20,seq_end)

   # dict to store best mapped read from each end
   lterminal_reads={}
   rterminal_reads={}

   for lread in left_reads:
      query_name = lread.query_name
      cigar = lread.cigarstring

      if lread.query_sequence is None: # skip empty reads
          continue
       
      # Check if the read is mapped multiple times and use
      # the read that maps to most bases
      if query_name in lterminal_reads:
         prior_read = lterminal_reads[query_name]
         prior_cigar = prior_read.cigarstring

         # Compare CIGAR strings to keep the one that maps more bases
         if cigar_maps_more_bases(cigar, prior_cigar):
            lterminal_reads[query_name] = lread
      else:
         lterminal_reads[query_name] = lread

   for rread in right_reads:
      query_name = rread.query_name
      cigar = rread.cigarstring

      if rread.query_sequence is None: # skip empty reads
         continue
      
      # Check if the read is mapped multiple times and use
      # the read that maps to most bases
      if query_name in rterminal_reads:
         prior_read = rterminal_reads[query_name]
         prior_cigar = prior_read.cigarstring

         # Compare CIGAR strings to keep the one that maps more bases
         if cigar_maps_more_bases(cigar, prior_cigar):
               rterminal_reads[query_name] = rread
      else:
         
         rterminal_reads[query_name] = rread

   # Write all fetched reads to a new file
   lterminal_file = pysam.AlignmentFile(loutput_handle,"w",template=input)
   for read in lterminal_reads.values():
      lterminal_file.write(read)
   lterminal_file.close()

   rterminal_file = pysam.AlignmentFile(routput_handle,"w",template=input)
   for read in rterminal_reads.values():
      rterminal_file.write(read)
   rterminal_file.close()

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

def is_map_empty(file_path):
    "Checks if a bam-file is empty by trying to fetch the first read"
    # Open the alignment file
    with pysam.AlignmentFile(file_path, "rb") as alignment_file:
        # Try to fetch the first read
        try:
            next(alignment_file)
            return False  # Alignment is not empty
        except StopIteration:
            return True  # Alignment is empty

def check_for_empty_cons(file_path):
   "Checks if a map was produced by an empty consensus"
   with pysam.AlignmentFile(file_path, "rb") as alignment_file:
      reads = list(alignment_file)  # Load all reads into a list
      
      # Check if there is exactly one read
      if len(reads) == 1:
         read = reads[0]
         # Check if the read is unmapped and has no sequence
         if read.is_unmapped and (not read.seq or read.seq == "*"):
            return True  # Only one unmapped read with no sequence
      return False  # Either more reads, or the read does not meet the conditions

# Example usage:

def stich_telo(ref,left_map,right_map,outfile,logout="consensus.log.txt"):
   """Writes a fasta-file extended with sequence based on a left- and right-side .bam-file. 
   Also Adds this information to a log file"""

   left_log_mes=""
   # Check if an empty left consensus was used to generate the map:
   if check_for_empty_cons(left_map):
      # Make an empty seq list to enable errors later on
      left_seqs=[]
      left_log_mes="#No consensus produced for left-side end. Likely, no reads extends the assembly. "
   else:
      # extract left cons-to-stich
      l_sam_in = pysam.AlignmentFile(left_map, "r")
      left_seqs =[]
      start_clip = r"^(\d+)S"
      # filter away mapping at right side
      cons_at_left = [read for read in l_sam_in if read.reference_start<1000]
      
      # Get the sequence extending beyond the genome
      for read in cons_at_left:
         lmatch = re.match(start_clip,read.cigarstring)
         if lmatch:
            clip_num=int(lmatch.group(1)) # digits are retrieve via .group
            seq = read.query_sequence[0:(clip_num-read.reference_start)] # Adjust for if more than just overhanging bases are soft-clipped
            left_seqs.append(seq)
      l_sam_in.close()
   right_log_mes=""
   # Check if an empty left consensus was used to generate the map:
   if check_for_empty_cons(right_map):
      right_seqs=[]
      right_log_mes="#No consensus produced for right-side end. Likely, no reads extends the assembly."
   else:
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
      logging.info("Left consensus does not extend genome")
   else:
      left_cons = SeqRecord(Seq(left_seqs[0]),id="left_cons")
      logging.info(f"Left consensus is {len(left_cons)}")
   if len(right_seqs)==0:
      right_cons=SeqRecord(Seq("")) # if it is empty make an empty seqrecord to avoid errors in joining later
      logging.info("Right cons does not extend genome")
   else:
      right_cons = SeqRecord(Seq(right_seqs[0]),id="right_cons")
   
      logging.info(f"Right consensus is {len(right_cons)}")
   new_genome = left_cons+genome+right_cons
   new_genome.id="Reference_with_consensus_attached"
   new_genome.description="" 
   SeqIO.write(new_genome,outfile,"fasta")
   SeqIO.write(right_cons,"tmp.right.fasta","fasta")
   SeqIO.write(left_cons,"tmp.left.fasta","fasta")
   
   # Create log of consensus length
   log = open(logout,"w")
   log.write("==============================================================================")
   log.write("\nINTIAL CONSENSUS")
   log.write("\n==============================================================================")
   log_content ="\nleft_cons:{}\tright_consensus:{}".format(len(left_cons),len(right_cons))
   comment_mes ="\n" + "\n".join([left_log_mes,right_log_mes]) 
   log_content=log_content + comment_mes
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

   # Get reference name from BAM file
   reference_name = bam_in.get_reference_name(0)
   coverage_count=bam_in.count_coverage(reference_name,start=position,stop=position+1,read_callback="nofilter",quality_threshold=1)
   A_num=coverage_count[0][0]
   C_num=coverage_count[1][0]
   G_num=coverage_count[2][0]
   T_num=coverage_count[3][0]
   cov=A_num+C_num+G_num+T_num
   
   if fasta_file.seq[position].upper()=="N":
      matching_bases=0
   elif fasta_file.seq[position].upper()=="A":
      matching_bases=A_num
   elif fasta_file.seq[position].upper()=="C":
      matching_bases=C_num
   elif fasta_file.seq[position].upper()=="G":
      matching_bases=G_num
   elif fasta_file.seq[position].upper()=="T":
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

         if cov>=cov_thres and (match/cov)>ratio_thres:
            index_start=pos

            break
      except TypeError: # if no reads are mapped
         continue
   
   #trim end/right
   for pos in range(fasta_end,fasta_end-right_len,-1):

      try:
         cov, match = get_support_info(sorted_bam_file,genome,pos,qual_thres)

         if cov>=cov_thres and (match/cov)>ratio_thres:
            index_end=pos

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
   log.write("\n==============================================================================")
   log.write("\nCONSENSUS TRIMMING")
   log.write("\n==============================================================================")
   log.write("\nRule: Trimmed until coverage>=5 and 70% read matching position:")
   log.write(log_message)
   log.close()
   SeqIO.write(trimmed_fasta,output_handle,"fasta")

def trim_by_map_illumina(genome, sorted_bam_file, output_handle,cons_log, cov_thres=1,ratio_thres=0.7, qual_thres=30):
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

         if cov>=cov_thres and (match/cov)>ratio_thres:
            index_start=pos

            break
      except TypeError: # if no reads are mapped
         continue
   
   #trim end/right
   for pos in range(fasta_end,fasta_end-right_len,-1):

      try:
         cov, match = get_support_info(sorted_bam_file,genome,pos,qual_thres)

         if cov>=cov_thres and (match/cov)>ratio_thres:
            index_end=pos

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
   log.write("\n==============================================================================")
   log.write("\nCONSENSUS TRIMMING")
   log.write("\n==============================================================================")
   log.write("\nRule: Trimmed until coverage Q_score>Q30")
   log.write(log_message)
   log.close()
   SeqIO.write(trimmed_fasta,output_handle,"fasta")

def generate_support_log(genome, qc_bam_file, output_handle):
   "Generates a coverage log for the ends of the genome"
   #trim start/left-side
   
   fasta = SeqIO.read(genome,"fasta")
   fasta_end = len(fasta.seq)-1 # subtract one to make it 0-indexed

   # Generate log of coverage at all positions
   with open(output_handle,"a") as log:

      for pos in range(0,fasta_end):
         
         try:
            cov, match = get_support_info(bam_file=qc_bam_file,
                                          genome=genome,
                                          position=pos,
                                          qual_threshold=1)
            
            print(pos,cov,match)
            log.write(pos,cov,match)
         except TypeError: # if no reads are mapped
            continue
   
   # Filter log for relevant positions

if __name__ == '__main__':
   print("testing module function")
   sam_to_readpair("test_files/left.sam",
                   "test_files/NBC_00001_FDMS210417786-1a_HMWN2DSX2_L3_1.fq.gz",
                   "test_files/NBC_00001_FDMS210417786-1a_HMWN2DSX2_L3_2.fq.gz",
                   "test_out1.fq",
                   "test_out2.fq")
   
   
