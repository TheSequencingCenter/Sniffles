# Author:  Richard Casey
# Date:    07-10-2024 (DD-MM-YYYY)
# Purpose: Utility function to delete all files in subdirectories.  Use with caution.

# Standard library imports
import os
import shutil

def clear_files():
    """
    Delete all files in subdirectories. 
    Use to set initial conditions before starting application.
    Use with caution.
    """

    # delete existing bam files
    bam_files_dir = os.environ.get('BAM_FILES_DIR')

    if bam_files_dir:
        for filename in os.listdir(bam_files_dir):
            file_path = os.path.join(bam_files_dir, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
            
    # delete existing fasta files
    fasta_files_dir = os.environ.get('FASTA_FILES_DIR')

    if fasta_files_dir:
        for filename in os.listdir(fasta_files_dir):
            file_path = os.path.join(fasta_files_dir, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)

    # delete existing fastq files
    fastq_files_dir = os.environ.get('FASTQ_FILES_DIR')

    if fastq_files_dir:
        for filename in os.listdir(fastq_files_dir):
            file_path = os.path.join(fastq_files_dir, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)

    # delete existing ref genome index files with suffix .fai and .gzi
    ref_files_dir = os.environ.get('REF_FILES_DIR')

    if ref_files_dir:
        for filename in os.listdir(ref_files_dir):
            if filename.endswith('.fai') or filename.endswith('.gzi'):
                file_path = os.path.join(ref_files_dir, filename)
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)

    # delete existing sniffles output files
    sniffles_files_dir = os.environ.get('SNIFFLES_FILES_DIR')

    if sniffles_files_dir:
        for filename in os.listdir(sniffles_files_dir):
            file_path = os.path.join(sniffles_files_dir, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)