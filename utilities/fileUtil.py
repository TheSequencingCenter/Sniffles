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

    # delete existing pod5 files
    pod5_files_dir = os.environ.get('POD5_FILES_DIR')

    if pod5_files_dir:
        for filename in os.listdir(pod5_files_dir):
            file_path = os.path.join(pod5_files_dir, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)

    # delete existing fast5 files
    fast5_files_dir = os.environ.get('FAST5_FILES_DIR')

    if fast5_files_dir:
        for filename in os.listdir(fast5_files_dir):
            file_path = os.path.join(fast5_files_dir, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)

    # delete existing bam files
    bam_files_dir = os.environ.get('BAM_FILES_DIR')

    if bam_files_dir:
        for filename in os.listdir(bam_files_dir):
            file_path = os.path.join(bam_files_dir, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)

    # delete existing sniffles output files
    sniffles_files_dir = os.environ.get('SNIFFLES_FILES_DIR')

    if sniffles_files_dir:
        for filename in os.listdir(sniffles_files_dir):
            file_path = os.path.join(sniffles_files_dir, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)