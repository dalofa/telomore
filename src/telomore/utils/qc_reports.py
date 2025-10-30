"""
Functions for generating useful QC metrics from the telomore script.
"""

import csv
import os
import tempfile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam

from .cmd_tools import map_and_sort, map_and_sort_illumina
from .fasta_tools import (
    cat_and_derep_fastq,
    check_fastq_order,
    dereplicate_fastq,
    merge_fasta,
)
from .map_tools import sam_to_fastq, sam_to_readpair


def qc_map(
    extended_assembly: str, left: str, right: str, output_handle: str, t=1
) -> None:
    """Collect terminal reads previously identified
    and maps them against the extended assembly"""
    # The file has to be mode=w to create the correct type of object
    # for sam_to_fastq
    with tempfile.NamedTemporaryFile(
        suffix='.fastq', delete=False, mode='w'
    ) as temp_fastq:
        temp_fastq_path = temp_fastq.name
        sam_to_fastq(left, temp_fastq)
        sam_to_fastq(right, temp_fastq)
    dereplicate_fastq(fastq_in=temp_fastq_path, fastq_out=temp_fastq_path)
    map_and_sort(extended_assembly, temp_fastq_path, output_handle, t)
    os.remove(temp_fastq_path)


def qc_map_illumina(
    extended_assembly: str,
    left_sam: str,
    right_sam: str,
    fastq_in1: str,
    fastq_in2: str,
    output_handle: str,
    t=1,
) -> None:
    """Collect terminal reads previously identified and maps them against extended assembly."""
    # get left paired read
    with (
        tempfile.TemporaryDirectory()
    ) as temp_dir:  # ensures files are deleted after usage
        # Create multiple temporary files in the temporary directory
        l_tmp1 = os.path.join(temp_dir, 'terminal_left_reads_1.fastq')
        l_tmp2 = os.path.join(temp_dir, 'terminal_left_reads_2.fastq')
        r_tmp1 = os.path.join(temp_dir, 'terminal_right_reads_1.fastq')
        r_tmp2 = os.path.join(temp_dir, 'terminal_right_reads_2.fastq')
        a_tmp1 = os.path.join(temp_dir, 'all_terminal_reads_1.fastq')
        a_tmp2 = os.path.join(temp_dir, 'all_terminal_reads_2.fastq')

        sam_to_readpair(
            sam_in=left_sam,
            fastq_in1=fastq_in1,
            fastq_in2=fastq_in2,
            fastq_out1=l_tmp1,
            fastq_out2=l_tmp2,
        )

        # get right paired read

        sam_to_readpair(
            sam_in=right_sam,
            fastq_in1=fastq_in1,
            fastq_in2=fastq_in2,
            fastq_out1=r_tmp1,
            fastq_out2=r_tmp2,
        )

        # collect the paired read files:
        cat_and_derep_fastq(fastq_in1=l_tmp1, fastq_in2=r_tmp1, fastq_out=a_tmp1)

        cat_and_derep_fastq(fastq_in1=l_tmp2, fastq_in2=r_tmp2, fastq_out=a_tmp2)

        if check_fastq_order(a_tmp1, a_tmp2):
            map_and_sort_illumina(
                reference=extended_assembly,
                read1=a_tmp1,
                read2=a_tmp2,
                output=output_handle,
                threads=t,
            )
        else:
            raise Exception('FASTQ files are not properly paired or ordered.')


def cons_genome_map(
    left_cons: str, right_cons: str, polished_genome: str, output_handle: str, t=1
) -> None:
    """Collect consensus and maps them against the polished genome"""
    merge_fasta(left_cons, right_cons, 'all_cons.fasta')
    map_and_sort(polished_genome, 'all_cons.fasta', output_handle, t)


def cons_cons_map(
    left_cons: str, right_cons: str, output_handle: str, t: int = 1
) -> None:
    """Collect consensus and maps them against each other"""
    map_and_sort(left_cons, right_cons, output_handle, t)


def cons_length(cons_file: str, output_handle: str, offset: int = 100) -> None:
    """Write the length of the consensus sequences to a tsv file."""
    cons_file = SeqIO.parse(cons_file, 'fasta')
    header = ['seq_id', 'end_cons', 'full_cons']
    tsv_log = []
    tsv_log.append(header)

    for record in cons_file:
        seq_id = record.id
        seq_len = len(record)
        gen_len = int(seq_len) - offset
        tsv_log.append([seq_id, gen_len, seq_len])

    with open(output_handle, 'w', newline='') as tsv_file:
        writer = csv.writer(tsv_file, delimiter='\t')
        writer.writerows(tsv_log)


def map_to_depth(bam_file: str, output_handle: str) -> None:
    """Map the depth of coverage to a bam file and write it to a file."""
    pysam.depth('-aa', bam_file, '-o', output_handle)


def finalize_log(log: str, right_fasta: str, left_fasta: str) -> None:
    """Function to finalize the log by prepending the final information about consensus added."""
    file = open(log)
    log_cont = file.readlines()
    file.close()
    # Org lengths of consensus added
    length_lines = log_cont[3]
    left_len = int(length_lines.split('\t')[0].split(':')[1])
    right_len = int(length_lines.split('\t')[1].split(':')[1])

    # get the number of bases trimmed off
    trim_left = log_cont[-2].split(' ')[-1]
    trim_right = log_cont[-1].split(' ')[-1]
    left_seq = SeqIO.read(left_fasta, 'fasta')
    right_seq = SeqIO.read(right_fasta, 'fasta')

    if trim_left.rstrip() == 'rejected':
        new_left = 'rejected'
        left_seq = SeqRecord(Seq(''))
    else:
        new_left = left_len - int(trim_left)
        left_seq = left_seq[int(trim_left) :]
    if trim_right.rstrip() == 'rejected':
        new_right = 'rejected'
        right_seq = SeqRecord(Seq(''))
    else:
        new_right = right_len - int(trim_right)
        right_seq = right_seq[0:new_right]

    final_lengths = 'left_cons:{}\tright_consensus:{}'.format(new_left, new_right)

    # write to log file
    file = open(log, 'w')
    file.write(
        '=============================================================================='
    )
    file.write('\nFINAL GENOME EXTENSION')
    file.write(
        '\n==============================================================================\n'
    )
    file.write(final_lengths)
    file.write('\n>left_cons\n')
    file.write(str(left_seq.seq))
    file.write('\n>right_cons\n')
    file.write(str(right_seq.seq))
    file.write('\n')
    for line in log_cont:
        file.write(line)
    file.write(
        '==============================================================================\n'
    )
    file.close()
