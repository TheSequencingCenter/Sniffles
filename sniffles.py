# Author:  Richard Casey
# Date:    07-24-2024 (DD-MM-YYYY)
# Purpose: Human genome structural variant caller.
# Notes:   In production mode, this application requires Nvidia GPU's (H100, A100, V100).

# Standard library imports
import argparse
import os
import subprocess
import sys

# Third-party imports
import boto3
from   botocore.exceptions import ClientError
from   botocore.exceptions import NoCredentialsError
from   botocore.exceptions import PartialCredentialsError

# Local application/library specific imports
from utilities.fileUtil   import clear_files
from utilities.loggerUtil import logger

def run_command(command: str) -> None:
    """Run a shell command and print the output."""
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        logger.info(f"{result.stdout}")
        logger.info(f"{result.stderr}")
    except subprocess.CalledProcessError as e:
        logger.error(f"ERROR: run_command failed with error: {e}")
        logger.error(f"STDOUT: {e.stdout}")
        logger.error(f"STDERR: {e.stderr}")
        sys.exit(1)

def set_environment_variables(sample_name: str) -> None:
    """Set environment variables."""
    try:
        # check which server is hosting the code 
        home_dir = os.path.expanduser("~")                             # get home directory
        if home_dir   == '/home/seqcenter':
            os.environ['BASE_DIR'] = '/home/seqcenter/Sniffles'        # for development ubuntu
        elif home_dir == '/home/ubuntu':
            os.environ['BASE_DIR'] = '/home/ubuntu/Sniffles'           # for production ubuntu on AWS
        else:
            logger.error("ERROR: unknown home directory setting")
            sys.exit(1)
    except ValueError as e:
        logger.error(f"Error: an unexpected error occured{e}")
        sys.exit(1)

    # MAIN DIRS
    os.environ['INPUT_DIR']          = f"{os.environ['BASE_DIR']}/data/INPUTS"
    os.environ['OUTPUT_DIR']         = f"{os.environ['BASE_DIR']}/data/OUTPUTS"

    # POD5
    os.environ['POD5_FILES_DIR']     = f"{os.environ['INPUT_DIR']}/POD5_FILES"

    # FAST5
    os.environ['FAST5_FILES_DIR']    = f"{os.environ['INPUT_DIR']}/FAST5_FILES"

    # # BAM
    os.environ['BAM_FILES_DIR']      = f"{os.environ['INPUT_DIR']}/BAM_FILES"
    os.environ['BAM_FILE']           = f"{os.environ['INPUT_DIR']}/BAM_FILES/{sample_name}.bam"
    os.environ['BAM_SORTED_FILE']    = f"{os.environ['INPUT_DIR']}/BAM_FILES/{sample_name}.sorted.bam"

    # REFERENCE 
    os.environ['REF_FILES_DIR']      = f"{os.environ['INPUT_DIR']}/REF_FILES"
    os.environ['REF_FILE']           = f"{os.environ['INPUT_DIR']}/REF_FILES/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz"

    # # MODELS
    os.environ['DORADO_MODELS']      = f"{os.environ['BASE_DIR']}/dorado_models"

    # SNIFFLES
    os.environ['SNIFFLES_FILES_DIR'] = f"{os.environ['OUTPUT_DIR']}/sniffles_output"
    os.environ['SNIFFLES_VCF_FILE']  = f"{os.environ['SNIFFLES_FILES_DIR']}/sniffles.vcf.gz"

    # THREADS
    os.environ['THREADS']            = '14'

def clear_all_files() -> None:
    """Clear all files in the specified directories."""
    logger.info("Clear files...")
    clear_files()

def check_s3_bucket_exists(bucket_name: str) -> bool:
    """Check if the specified S3 bucket exists."""
    s3 = boto3.client('s3')
    try:
        s3.head_bucket(Bucket=bucket_name)
        return True
    except ClientError as e:
        if e.response['Error']['Code'] == '404':
            return False
        else:
            raise
    except NoCredentialsError:
        logger.error("ERROR: AWS credentials not found.")
        sys.exit(1)
    except PartialCredentialsError:
        logger.error("ERROR: Incomplete AWS credentials.")
        sys.exit(1)

def copy_sample_files_S3_to_EC2(bucket_name: str, sample_name: str) -> None:
    """Copy sample files from S3 to EC2.
    
     Parameters:
       --recursive : Copy files recursively.  This is a required parameter.
    """
    command = (
        # f"aws s3 cp s3://{bucket_name}/{sample_name}/ "
        f"aws s3 cp s3://{bucket_name}/ "
        f"{os.environ['POD5_FILES_DIR']} "
        "--recursive"
    )
    run_command(command)

def convert_fast5_to_pod5() -> None:
    """Convert FAST5 file to POD5 file.

    Parameters:
       -o : POD5 output directory.
       -t : Number of threads
       -f : force overwrite existing POD5 files 
    """
    command = (
        "pod5 convert fast5 "
        f"-o {os.environ['POD5_FILES_DIR']} "
        f"-t {os.environ['THREADS']} "
        "-f "
        f"{os.environ['FAST5_FILES_DIR']}"
    )
    run_command(command)

def convert_pod5_to_bam() -> None:
    """Convert POD5 file to BAM file using dorado basecaller.
    Note: In production mode, set "-x" to "cuda:all" so it uses Nvidia GPU's (H100, A100, V100).
    
    Parameters:
        -x          : use cpu only or use gpu's. device string in format "cuda:0,...,N", "cuda:all", "metal", "cpu" etc.. [default: "cuda:all"]
        --batchsize : if 0, an optimal batchsize will be selected. batchsizes are rounded to the closest multiple of 64. This may only affect the amount of GPU RAM required to run. [default: 0]
        -v          : verbose mode 
    """
    command = (
        "dorado basecaller "
        "-x cpu "      
        # "-x cuda:all "
        # "--batchsize 64 "
        # "-v "
        f"--reference {os.environ['REF_FILE']} "
        f"            {os.environ['DORADO_MODELS']}/dna_r10.4.1_e8.2_400bps_fast@v4.2.0 "
        f"            {os.environ['POD5_FILES_DIR']} > "  # this parameter is either a single file name or a directory name (confusing)
        f"            {os.environ['BAM_FILE']}"
    )
    run_command(command)

def create_bam_index_file() -> None:
    """Create bam index file."""
    command = (
        "samtools index "
        f"-@ {os.environ['THREADS']} "
        f"   {os.environ['BAM_FILE']}"
    )
    run_command(command)

def sort_bam_file() -> None:
    """Sort bam file."""
    command = (
        "samtools sort "
        f"--threads {os.environ['THREADS']} "
        f"          {os.environ['BAM_FILE']} "
        f"-o        {os.environ['BAM_SORTED_FILE']}"
    )
    run_command(command)

def create_sorted_bam_index_file() -> None:
    """Create sorted bam index file."""
    command = (
        "samtools index "
        f"-@ {os.environ['THREADS']} "
        f"   {os.environ['BAM_SORTED_FILE']}"
    )
    run_command(command)

def run_sniffles() -> None:
    """
    Perform sturctural variant calling.

    Parameters:
        -i                : Input sorted BAM file.
        -v                : Output sniffles VCF file.
        -t                : Number of threads.
        --allow-overwrite : Allow overwrite of existing sniffles VCF file.
    """

    command = (
        "sniffles "
        f"-i {os.environ['BAM_SORTED_FILE']} "
        f"-v {os.environ['SNIFFLES_VCF_FILE']} "
        f"-t {os.environ['THREADS']} "
        "--allow-overwrite"
    )
    run_command(command)

# main entry point
if __name__ == "__main__":

    try:
        # set S3 bucket name and sample name
        parser = argparse.ArgumentParser()
        parser.add_argument("-b", "--bucketname", required=True, help="Enter a bucket name for this run.")
        parser.add_argument("-s", "--samplename", required=True, help="Enter a sample name for this run.")
        args = parser.parse_args()
            
        logger.info("Start Sniffles...")

        # set env vars
        try:
            logger.info("Set environment variables...")
            set_environment_variables(args.samplename)
        except Exception as e:
            logger.error(f"ERROR: Could not set environment variables {e}")
            sys.exit(1)

        # Check if S3 bucket exists
        try:
            if not check_s3_bucket_exists(args.bucketname):
                logger.error(f"ERROR: AWS S3 bucket '{args.bucketname}' does not exist.")
                sys.exit(1)
        except NoCredentialsError:
            logger.error("ERROR: AWS credentials not found.")
            sys.exit(1)
        except PartialCredentialsError:
            logger.error("ERROR: Incomplete AWS credentials.")
            sys.exit(1)
        except Exception as e:
            logger.error(f"ERROR: An error occurred while checking the AWS S3 bucket: {e}")
            sys.exit(1)

        # delete all files in subdirectories to set initial state for application
        try:
            logger.info("Clear files...")
            clear_files()
        except Exception as e:
            logger.error(f"ERROR: Failed to clear files: {e}")
            sys.exit(1)

        # 1. copy sample files from S3 seqcenter-samples bucket to EC2 POD5 directory
        try:
            logger.info("Copy sample POD5 files from S3 to EC2...")
            copy_sample_files_S3_to_EC2(args.bucketname, args.samplename)
        except Exception as e:
            logger.error(f"ERROR: Failed to copy sample POD5 files from S3 to EC2: {e}")
            sys.exit(1)

        # convert fast5 file to pod5 file
        # try:
        #     logger.info("Convert fast5 to pod5...")
        #     convert_fast5_to_pod5()
        # except Exception as e:
        #     logger.error(f"ERROR: Failed to convert fast5 to pod5: {e}")
        #     sys.exit(1)

        # 2. convert pod5 file to bam file
        try:
            logger.info("Convert pod5 to bam...")
            convert_pod5_to_bam()
        except Exception as e:
            logger.error(f"ERROR: Failed to convert pod5 to bam: {e}")
            sys.exit(1)

        # 3. sort bam file
        try:
            logger.info("Sort bam file...")
            sort_bam_file()
        except Exception as e:
            logger.error(f"ERROR: Failed to sort bam file: {e}")
            sys.exit(1)

        # 4. create sorted bam index file
        try:
            logger.info("Create sorted bam index file...")
            create_sorted_bam_index_file()
        except Exception as e:
            logger.error(f"ERROR: Failed to create sorted bam index file: {e}")
            sys.exit(1)

        # 5. perform structural variant calling
        try:
            logger.info("Perform structural variant calling...")
            run_sniffles()
        except Exception as e:
            logger.error(f"ERROR: Failed to perform structural variant calling: {e}")
            sys.exit(1)

        logger.info("Finish Sniffles...")

    except Exception as e:
        logger.error(f"ERROR: an error occurred in main: {e}")
        sys.exit(1)