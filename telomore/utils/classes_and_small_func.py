"""
Class for handling files related to each replicon in a file
"""
import os
import shutil

# Class
# Replicon class
class Replicon:
    """Class to handle files related to each replicon in a file"""
    def __init__(self, name: str, org_fasta):
        self.name = name
        self.org_fasta = org_fasta

    # Map files
        self.org_map = f"{self.name}_map.bam"
        self.org_map_index = f"{self.name}_map.bam.bai"

        # Filtered files
        self.left_sam = f"{self.name}_left.sam"
        self.left_filt = f"{self.name}_left_filtered"
        self.left_filt_sam = f"{self.name}_left_filtered.sam"
        self.left_filt_fq = f"{self.name}_left_filtered.fastq"

        self.right_sam = f"{self.name}_right.sam"
        self.right_filt = f"{self.name}_right_filtered"
        self.right_filt_sam = f"{self.name}_right_filtered.sam"
        self.right_filt_fq = f"{self.name}_right_filtered.fastq"

    # Consensus files
        # left
        self.l_cons_out =f"rev_{self.name}_left_cons.fasta"
        self.l_cons_final_out = f"{self.name}_left_cons.fasta"
        self.l_cons_alignment=f"{self.l_cons_out}.aln"
        self.revcomp_out = f"rev_{self.left_filt_fq}"
        # right
        self.r_cons_final_out = f"{self.name}_right_cons.fasta"
        self.r_cons_alignment=f"{self.r_cons_final_out}.aln"

    # Extension files
        self.contig_fasta = f"{name}.fasta"

        self.cons_log_np_out= f"{self.name}_telomore_ext_np.log"
        self.cons_log_ill_out = f"{self.name}_telomore_ill_ext.log"

        # Truncated contig which discard alternative mapping points
        self.trunc_left_fasta = f"{self.name}_trunc_left.fa"
        self.trunc_right_fasta = f"{self.name}_trunc_right.fa"

        # Maps on trunc fasta
        self.l_map_out = f"{self.name}_left_map.bam"
        self.r_map_out = f"{self.name}_right_map.bam"
        self.l_map_out_index = f"{self.l_map_out}.bai"
        self.r_map_out_index = f"{self.r_map_out}.bai"

        # Extended assembly
        self.stitch_out = f"{self.name}_telomore_untrimmed.fasta"
        self.stitch_left_fasta = f"{self.name}_left.fasta"
        self.stitch_right_fasta = f"{self.name}_right.fasta"

    # Trim files
        self.trim_map = f"{self.name}_telomore_untrimmed.bam"
        self.trim_map_index = f"{self.trim_map}.bai"
        self.trim_out = f"{self.name}_telomore_extended.fasta"

    # QC_files
        self.qc_out = f"{self.name}_telomore_QC.bam"
        self.qc_out_index = f"{self.qc_out}.bai"

    def cleanup_tmp_files(self):
        """Clean up temporary files"""
        tmp_files=[self.org_map,
                   self.org_map_index,
                   self.left_sam,
                   self.left_filt_sam,
                   self.left_filt_fq,
                   self.right_sam,
                   self.right_filt_sam,
                   self.right_filt_fq,
                   self.l_cons_out,
                   self.l_cons_final_out,
                   self.l_cons_alignment,
                   self.revcomp_out,
                   self.r_cons_final_out,
                   self.r_cons_alignment,
                   self.contig_fasta,
                   self.trunc_left_fasta,
                   self.trunc_right_fasta,
                   self.l_map_out,
                   self.l_map_out_index,
                   self.r_map_out_index,
                   self.r_map_out,
                   self.stitch_left_fasta,
                   self.stitch_right_fasta,
                   self.trim_map,
                   self.trim_map_index
        ]
        for path in tmp_files:
            if os.path.exists(path):
                os.remove(path)

    def mv_files(self,folder,mode):
        """Move files that should be kept as part of the output"""
        keep_files= [self.stitch_out,
                     self.trim_out,
                     self.qc_out,
                     self.qc_out_index]

        for file in keep_files:
            shutil.move(src=file,
            dst = os.path.join(folder,file))
        if mode=="nanopore":
            shutil.move(src=self.cons_log_np_out,
                        dst=os.path.join(folder,self.cons_log_np_out))
        elif mode=="illumina":
            shutil.move(src=self.cons_log_ill_out,
            dst=os.path.join(folder,self.cons_log_ill_out))
