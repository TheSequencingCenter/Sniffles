# Author:  Richard Casey
# Date:    07-24-2024 (DD-MM-YYYY)
# Purpose: Human genome structural variant caller.

# Standard library imports
import argparse
import json
import os
import shutil
import subprocess
import sys

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

def set_environment_variables() -> None:
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

    # SAMPLES
    # sample_file_name = "test"
    # sample_file_name = "HG002_guppy422_2_GRCh38_no_alt"
    sample_file_name = "giab_2023_05"

    # MAIN DIRS
    os.environ['INPUT_DIR']              = f"{os.environ['BASE_DIR']}/data/INPUTS"
    os.environ['OUTPUT_DIR']             = f"{os.environ['BASE_DIR']}/data/OUTPUTS"

    # POD5
    os.environ['POD5_FILES_DIR']         = f"{os.environ['INPUT_DIR']}/POD5_FILES"
    os.environ['POD5_FILE']              = f"{os.environ['INPUT_DIR']}/POD5_FILES/{sample_file_name}.pod5"

    # FAST5
    os.environ['FAST5_FILES_DIR']        = f"{os.environ['INPUT_DIR']}/FAST5_FILES"
    os.environ['FAST5_FILE']             = f"{os.environ['INPUT_DIR']}/FAST5_FILES/{sample_file_name}.fast5"

    # FASTQ
    os.environ['FASTQ_FILES_DIR']        = f"{os.environ['INPUT_DIR']}/FASTQ_FILES"
    os.environ['FASTQ_FILE']             = f"{os.environ['INPUT_DIR']}/FASTQ_FILES/{sample_file_name}.fastq"

    # FASTA
    os.environ['FASTA_FILES_DIR']        = f"{os.environ['INPUT_DIR']}/FASTA_FILES"
    os.environ['FASTA_FILE']             = f"{os.environ['INPUT_DIR']}/FASTA_FILES/{sample_file_name}.fasta"

    # # BAM
    os.environ['BAM_FILES_DIR']          = f"{os.environ['INPUT_DIR']}/BAM_FILES"
    os.environ['BAM_FILE']               = f"{os.environ['INPUT_DIR']}/BAM_FILES/{sample_file_name}.bam"
    os.environ['BAM_SORTED_FILE']        = f"{os.environ['INPUT_DIR']}/BAM_FILES/{sample_file_name}.sorted.bam"

    # REFERENCE 
    os.environ['REF_FILES_DIR']          = f"{os.environ['INPUT_DIR']}/REF_FILES"
    os.environ['REF_FILE']               = f"{os.environ['INPUT_DIR']}/REF_FILES/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz"

    # # MODELS
    os.environ['DORADO_MODELS']          = f"{os.environ['BASE_DIR']}/dorado_models"

    # SNIFFLES
    os.environ['SNIFFLES_FILES_DIR']    = f"{os.environ['OUTPUT_DIR']}/sniffles_output"
    os.environ['SNIFFLES_VCF_FILE']      = f"{os.environ['SNIFFLES_FILES_DIR']}/sniffles.vcf.gz"

    # THREADS
    os.environ['THREADS']                = '14'

def load_mappings(environment):
    """Manage input-output files. Load the mappings for the specified environment from the JSON file."""
    mappings_file = f"data/io-mappings/{environment}_io_mappings.json"
    if os.path.exists(mappings_file):
        with open(mappings_file, 'r') as file:
            return json.load(file)
    else:
        return {"io-mappings": []}

def save_mappings(environment, mappings):
    """Manage input-output files. Save the updated mappings to the JSON file for the specified environment."""
    mappings_file = f"data/io-mappings/{environment}_io_mappings.json"
    with open(mappings_file, 'w') as file:
        json.dump(mappings, file, indent=4)

def add_mapping(environment, input_file, output_file):
    """Manage input-output files. Add a new mapping of input and output files for the specified environment."""
    mappings = load_mappings(environment)
    mappings["io-mappings"].append({"input": input_file, "output": output_file})
    save_mappings(environment, mappings)

def get_output_file(environment, input_file):
    """Manage input-output files. Retrieve the output file corresponding to the given input file for the specified environment."""
    mappings = load_mappings(environment)
    for mapping in mappings["io-mappings"]:
        if mapping["input"] == input_file:
            return mapping["output"]
    return None

def copy_sample_files_S3_to_EC2() -> None:
    """Copy sample files from S3 to EC2.
    
    Parameters:
        --recursive : Copy all files recursively.
    """
    command = (
        "aws s3 cp s3://seqcenter-samples/samples/ "
        f"{os.environ['POD5_FILES_DIR']} "
        "--recursive"
    )
    run_command(command)

def convert_pod5_to_bam() -> None:
    """Convert POD5 file to BAM file using dorado basecaller.
    
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

def convert_fast5_to_pod5() -> None:
    """Convert FAST5 file to POD5 file using dorado basecaller.

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

def convert_fast5_to_bam() -> None:
    """Convert FAST5 file to BAM file using dorado basecaller.

    Parameters:time docker run \
    -v $(pwd):/data \
    pepper_deepvariant \
    margin phase \
    /data/HG002_guppy422_2_GRCh38_no_alt.bam \
    /data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz \
    /data/PEPPER_VARIANT_FULL.vcf.gz \
    /opt/margin_dir/params/phase/allParams.phase_vcf.ont.json \
    -o /data/margin_output/HG002_guppy422_2_GRCh38_no_alt \
    -t 14

        -x : cpu only
    """
    command = (
        "dorado basecaller "
        "-x cpu "
        f"{os.environ['DORADO_MODELS']}/dna_r10.4.1_e8.2_400bps_fast@v4.2.0 "
        f"{os.environ['FAST5_FILE']} > " 
        f"{os.environ['BAM_FILE']}"
    )
    run_command(command)

def convert_pod5_to_fastq() -> None:
    """Convert POD5 file to FASTQ file using dorado basecaller.

    Parameters:
        -x : cpu only
        --emit-fastq
    """

    # delete existing fastq files
    fastq_file_dir = os.environ.get('FASTQ_FILES_DIR')

    if fastq_file_dir:
        for filename in os.listdir(fastq_file_dir):
            file_path = os.path.join(fastq_file_dir, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)

    command = (
        "dorado basecaller "
        "-x cpu "
        "--emit-fastq "
        f"{os.environ['DORADO_MODELS']}/dna_r10.4.1_e8.2_400bps_fast@v4.2.0 "
        f"{os.environ['POD5_FILE']} > " 
        f"{os.environ['FASTQ_FILE']}"
    )
    run_command(command)

def convert_fastq_to_fasta() -> None:
    """Convert FASTQ file to FASTA file"""

    # delete existing fasta files
    fasta_file_dir = os.environ.get('FASTA_FILES_DIR')

    if fasta_file_dir:
        for filename in os.listdir(fasta_file_dir):
            file_path = os.path.join(fasta_file_dir, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)

    command = (
        "seqtk seq "
        f"-a {os.environ['FASTQ_FILE']} > " 
        f"   {os.environ['FASTA_FILE']}"
    )
    run_command(command)

def create_fasta_index_file() -> None:
    """Create FASTA index file."""
    command = (
        "samtools faidx "
        f"{os.environ['FASTA_FILE']}"
    )
    run_command(command)

def create_ref_genome_index_file() -> None:
    """Create reference genome index file.
    Reference genome file must be compressed with bgzip, not gzip."""

    # delete existing ref genome index files with suffix .fai and .gzi
    ref_files_dir = os.environ.get('REF_FILES_DIR')

    if ref_files_dir:
        for filename in os.listdir(ref_files_dir):
            if filename.endswith('.fai') or filename.endswith('.gzi'):
                file_path = os.path.join(ref_files_dir, filename)
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)

    command = (
        "samtools faidx "
        f"{os.environ['REF_FILE']}"
    )
    run_command(command)

def align_fasta_to_reference() -> None:
    """Align FASTA file to reference genome using dorado basecaller."""

    # delete existing bam files
    bam_files_dir = os.environ.get('BAM_FILES_DIR')

    if bam_files_dir:
        for filename in os.listdir(bam_files_dir):
            file_path = os.path.join(bam_files_dir, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)

    command = (
        "dorado aligner "
        f"-t {os.environ['THREADS']} "
        f"   {os.environ['REF_FILE']} "
        f"   {os.environ['FASTA_FILE']} > "
        f"   {os.environ['BAM_FILE']}" 
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

    # clear files before starting application
    parser = argparse.ArgumentParser(description="Human genome structural variant caller.")
    parser.add_argument("--clearfiles", action="store_true", help="Clear files before starting application.")
    args   = parser.parse_args()

    try:
        logger.info("Start Sniffles...")

        # set env vars
        logger.info("Set environment variables...")
        set_environment_variables()

        # delete all files in subdirectories if --clearfiles is specified. Use with caution.
        if args.clearfiles:
            logger.info("Clear files...")
            clear_files()

        # # Add new mappings
        # add_mapping("test",       "data/input/test/input3.txt",       "data/output/test/output3.txt")
        # add_mapping("production", "data/input/production/input3.txt", "data/output/production/output3.txt")
        
        # # Get output file for a given input file
        # output_file = get_output_file("test", "data/input/test/{sample_file_name}.pod5")
        # print(f"Output file for 'data/output/test/{sample_file_name}.pod5': {output_file}")
        # output_file = get_output_file("production", "data/input/production/input3.txt")
        # print(f"Output file for 'data/input/production/input3.txt': {output_file}")

        # 1. copy sample files from S3 seqcenter-samples bucket to EC2 POD5 directory
        try:
            logger.info("Copy sample files from S3 to EC2...")
            copy_sample_files_S3_to_EC2()
        except Exception as e:
            logger.error(f"ERROR: Failed to copy files from S3 to EC2: {e}")

        # 2. convert pod5 file to bam file
        try:
            logger.info("Convert pod5 to bam...")
            convert_pod5_to_bam()
        except Exception as e:
            logger.error(f"ERROR: Failed to convert pod5 to bam: {e}")

        # convert fast5 file to bam file
        # try:
        #     logger.info("Convert fast5 to bam...")
        #     convert_fast5_to_bam()
        # except Exception as e:
        #     logger.error(f"ERROR: Failed to convert fast5 to bam: {e}")

        # convert fast5 file to pod5 file
        # try:
        #     logger.info("Convert fast5 to pod5...")
        #     convert_fast5_to_pod5()
        # except Exception as e:
        #     logger.error(f"ERROR: Failed to convert fast5 to pod5: {e}")

        # convert pod5 file to fastq file
        # try:
        #     logger.info("Convert pod5 to fastq...")
        #     convert_pod5_to_fastq()
        # except Exception as e:
        #     logger.error(f"ERROR: Failed to convert pod5 to fastq: {e}")

        # convert fastq file to fasta file
        # try:
        #     logger.info("Convert fastq to fasta...")
        #     convert_fastq_to_fasta()
        # except Exception as e:
        #     logger.error(f"ERROR: Failed to convert fastq to fasta: {e}")

        # create fasta index file
        # try:
        #     logger.info("Create fasta index file...")
        #     create_fasta_index_file()
        # except Exception as e:
        #     logger.error(f"ERROR: Failed to create fasta index file: {e}")

        # create reference genome index file
        # try:
        #     logger.info("Create reference genome index file...")
        #     create_ref_genome_index_file()
        # except Exception as e:
        #     logger.error(f"ERROR: Failed to create reference genome index file: {e}")

        # align fasta file to reference genome
        # try:
        #     logger.info("Align fasta to reference genome...")
        #     align_fasta_to_reference()
        # except Exception as e:
        #     logger.error(f"ERROR: Failed to align fasta to reference genome: {e}")

        # 3. sort bam file
        try:
            logger.info("Sort bam file...")
            sort_bam_file()
        except Exception as e:
            logger.error(f"ERROR: Failed to sort bam file: {e}")

        # 4. create sorted bam index file
        try:
            logger.info("Create sorted bam index file...")
            create_sorted_bam_index_file()
        except Exception as e:
            logger.error(f"ERROR: Failed to create sorted bam index file: {e}")

        # 5. perform structural variant calling
        try:
            logger.info("Perform structural variant calling...")
            run_sniffles()
        except Exception as e:
            logger.error(f"ERROR: Failed to perform structural variant calling: {e}")

        logger.info("Finish Sniffles...")

    except Exception as e:
        logger.error(f"ERROR: an error occurred in main: {e}")
        sys.exit(1)
