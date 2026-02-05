"""Class for handling files related to each replicon."""

import os
import shutil


# Class
# Replicon class
class Replicon:
    """
    Manage file paths and operations for individual replicon extension.

    A Replicon object tracks all input, intermediate, and output files
    associated with extending a single linear contig (replicon). It provides
    organized file management and cleanup methods for the Telomore workflow.

    Parameters
    ----------
    name : str
        Identifier for the replicon (typically contig name from FASTA)
    org_fasta : str
        Path to the original multi-FASTA file containing this replicon

    Attributes
    ----------
    name : str
        Replicon identifier
    org_fasta : str
        Original FASTA file path

    Mapping files:
        org_map : str
            BAM file of reads mapped to original contig
        org_map_index : str
            BAI index for org_map

    Terminal read files:
        left_sam, right_sam : str
            SAM files with left/right terminal read alignments
        left_filt, right_filt : str
            Base paths for filtered reads
        left_filt_sam, right_filt_sam : str
            Filtered SAM files
        left_filt_fq, right_filt_fq : str
            Filtered FASTQ files

    Consensus files:
        l_cons_out : str
            Left consensus (reverse complement, temporary)
        l_cons_final_out : str
            Left consensus (final orientation)
        l_cons_alignment : str
            Alignment file for left consensus
        revcomp_out : str
            Reverse complement of left filtered reads
        r_cons_final_out : str
            Right consensus (final)
        r_cons_alignment : str
            Alignment file for right consensus

    Extension files:
        contig_fasta : str
            Extracted single-contig FASTA
        cons_log_np_out : str
            Extension log for Nanopore mode
        cons_log_ill_out : str
            Extension log for Illumina mode
        trunc_left_fasta, trunc_right_fasta : str
            Truncated contigs to prevent alternative mappings
        l_map_out, r_map_out : str
            BAM files of consensus mapped to truncated contigs
        l_map_out_index, r_map_out_index : str
            BAI indices for consensus maps
        stitch_out : str
            Extended assembly before trimming
        stitch_left_fasta, stitch_right_fasta : str
            Extracted consensus sequences
        trim_map : str
            BAM of QC reads mapped to untrimmed assembly
        trim_map_index : str
            BAI index for trim_map
        trim_out : str
            Final trimmed extended assembly

    QC files:
        qc_out : str
            Final QC BAM file
        qc_out_index : str
            BAI index for QC BAM

    Notes
    -----
    All file paths are automatically generated from the replicon name
    following a consistent naming convention. This ensures files are
    traceable and organized.

    The class provides methods for:
    - cleanup_tmp_files(): Remove intermediate processing files
    - mv_files(): Move final output files to designated directory

    Files are categorized as:
    - Temporary: Deleted after successful extension
    - Output: Moved to results directory for user
    """

    def __init__(self, name: str, org_fasta: str):
        self.name = name
        self.org_fasta = org_fasta

        # Map files
        self.org_map = f'{self.name}_map.bam'
        self.org_map_index = f'{self.name}_map.bam.bai'

        # Filtered files
        self.left_sam = f'{self.name}_left.sam'
        self.left_filt = f'{self.name}_left_filtered'
        self.left_filt_sam = f'{self.name}_left_filtered.sam'
        self.left_filt_fq = f'{self.name}_left_filtered.fastq'

        self.right_sam = f'{self.name}_right.sam'
        self.right_filt = f'{self.name}_right_filtered'
        self.right_filt_sam = f'{self.name}_right_filtered.sam'
        self.right_filt_fq = f'{self.name}_right_filtered.fastq'

        # Consensus files
        # left
        self.l_cons_out = f'rev_{self.name}_left_cons.fasta'
        self.l_cons_final_out = f'{self.name}_left_cons.fasta'
        self.l_cons_alignment = f'{self.l_cons_out}.aln'
        self.revcomp_out = f'rev_{self.left_filt_fq}'
        # right
        self.r_cons_final_out = f'{self.name}_right_cons.fasta'
        self.r_cons_alignment = f'{self.r_cons_final_out}.aln'

        # Extension files
        self.contig_fasta = f'{name}.fasta'

        self.cons_log_np_out = f'{self.name}_telomore_ext_np.log'
        self.cons_log_ill_out = f'{self.name}_telomore_ill_ext.log'

        # Truncated contig which discard alternative mapping points
        self.trunc_left_fasta = f'{self.name}_trunc_left.fa'
        self.trunc_right_fasta = f'{self.name}_trunc_right.fa'

        # Maps on trunc fasta
        self.l_map_out = f'{self.name}_left_map.bam'
        self.r_map_out = f'{self.name}_right_map.bam'
        self.l_map_out_index = f'{self.l_map_out}.bai'
        self.r_map_out_index = f'{self.r_map_out}.bai'

        # Extended assembly
        self.stitch_out = f'{self.name}_telomore_untrimmed.fasta'
        self.stitch_left_fasta = f'{self.name}_left.fasta'
        self.stitch_right_fasta = f'{self.name}_right.fasta'

        # Trim files
        self.trim_map = f'{self.name}_telomore_untrimmed.bam'
        self.trim_map_index = f'{self.trim_map}.bai'
        self.trim_out = f'{self.name}_telomore_extended.fasta'

        # QC_files
        self.qc_out = f'{self.name}_telomore_QC.bam'
        self.qc_out_index = f'{self.qc_out}.bai'

    def cleanup_tmp_files(self) -> None:
        """
        Remove temporary intermediate files after successful extension.

        Deletes all intermediate files that are not needed in the final output,
        including mapping files, filtered reads, consensus intermediates, and
        truncated references. Preserves only the final extended assemblies and
        QC files.

        Returns
        -------
        None
            Removes files from the filesystem

        Notes
        -----
        Files removed include:
        - Original mapping: org_map, org_map_index
        - Terminal read SAMs: left_sam, right_sam
        - Filtered reads: left_filt_sam, left_filt_fq, right_filt_sam, right_filt_fq
        - Consensus intermediates: l_cons_out, l_cons_final_out, l_cons_alignment,
          revcomp_out, r_cons_final_out, r_cons_alignment
        - Extracted/truncated contigs: contig_fasta, trunc_left_fasta, trunc_right_fasta
        - Consensus mappings: l_map_out, r_map_out, and their indices
        - Stitching intermediates: stitch_left_fasta, stitch_right_fasta
        - Trimming map: trim_map, trim_map_index

        Files preserved (not deleted):
        - stitch_out: Untrimmed extended assembly
        - trim_out: Final trimmed extended assembly
        - qc_out, qc_out_index: QC alignment files
        - cons_log_np_out or cons_log_ill_out: Extension logs

        Only deletes files that exist - missing files are silently skipped.
        Call this method after successful completion of extension workflow
        to reduce disk space usage.
        """
        tmp_files = [
            self.org_map,
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
            self.trim_map_index,
        ]
        for path in tmp_files:
            if os.path.exists(path):
                os.remove(path)

    def mv_files(self, folder: str, mode: str) -> None:
        """
        Move final output files to designated output directory.

        Relocates the essential output files (extended assemblies, QC BAM, and
        extension log) from the working directory to the specified output folder.
        The log file moved depends on the sequencing mode.

        Parameters
        ----------
        folder : str
            Path to destination directory for output files
        mode : str
            Sequencing technology mode: 'nanopore' or 'illumina'

        Returns
        -------
        None
            Moves files to the destination folder

        Raises
        ------
        FileNotFoundError
            If any of the required output files don't exist (implicitly from shutil.move)

        Notes
        -----
        Files moved for all modes:
        - stitch_out: Untrimmed extended assembly
        - trim_out: Final trimmed extended assembly
        - qc_out: QC alignment BAM
        - qc_out_index: QC alignment BAM index

        Mode-specific files:
        - If mode='nanopore': moves cons_log_np_out
        - If mode='illumina': moves cons_log_ill_out

        The destination folder must already exist. Files retain their
        original names in the destination directory.

        This method should be called after cleanup_tmp_files() to organize
        the final results while removing intermediate files from the working
        directory.
        """
        keep_files = [self.stitch_out, self.trim_out, self.qc_out, self.qc_out_index]

        for file in keep_files:
            shutil.move(src=file, dst=os.path.join(folder, file))
        if mode == 'nanopore':
            shutil.move(
                src=self.cons_log_np_out, dst=os.path.join(folder, self.cons_log_np_out)
            )
        elif mode == 'illumina':
            shutil.move(
                src=self.cons_log_ill_out,
                dst=os.path.join(folder, self.cons_log_ill_out),
            )
