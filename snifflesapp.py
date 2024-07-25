# Author:  Richard Casey
# Date:    07-24-2024 (DD-MM-YYYY)
# Purpose: Human genome structural variant caller.

# Standard library imports
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
        # logger.error(f"ERROR: run_command failed with error: {e}")
        # sys.exit(1)
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
        logger.error(f"Error: an unexpedted error occured{e}")
        sys.exit(1)

    # SAMPLE
    sample_file_name = "test"
    # sample_file_name = "HG002_guppy422_2_GRCh38_no_alt"

    # MAIN DIRS
    # os.environ['TRAIN_DIR']              = f"{os.environ['BASE_DIR']}/data/PEPPER_TRAINING"
    os.environ['INPUT_DIR']              = f"{os.environ['BASE_DIR']}/data/INPUTS"
    os.environ['OUTPUT_DIR']             = f"{os.environ['BASE_DIR']}/data/OUTPUTS"

    # POD5
    os.environ['POD5_FILE']              = f"{os.environ['INPUT_DIR']}/POD5_FILES/{sample_file_name}.pod5"

    # FAST5
    os.environ['FAST5_FILE']             = f"{os.environ['INPUT_DIR']}/FAST5_FILES/{sample_file_name}.fast5"

    # # FASTQ
    # os.environ['FASTQ_FILES_DIR']        = f"{os.environ['INPUT_DIR']}/FASTQ_FILES"
    # os.environ['FASTQ_FILE']             = f"{os.environ['INPUT_DIR']}/FASTQ_FILES/{sample_file_name}.fastq"

    # # FASTA
    # os.environ['FASTA_FILES_DIR']        = f"{os.environ['INPUT_DIR']}/FASTA_FILES"
    # os.environ['FASTA_FILE']             = f"{os.environ['INPUT_DIR']}/FASTA_FILES/{sample_file_name}.fasta"

    # # BAM
    os.environ['BAM_FILES_DIR']          = f"{os.environ['INPUT_DIR']}/BAM_FILES"
    os.environ['BAM_FILE']               = f"{os.environ['INPUT_DIR']}/BAM_FILES/{sample_file_name}.bam"
    os.environ['BAM_SORTED_FILE']        = f"{os.environ['INPUT_DIR']}/BAM_FILES/{sample_file_name}.sorted.bam"

    # # TRUTH VCF
    # os.environ['TRUTH_VCF_FILES_DIR']    = f"{os.environ['INPUT_DIR']}/TRUTH_VCF_FILES"
    # os.environ['TRUTH_VCF_FILE']         = f"{os.environ['INPUT_DIR']}/TRUTH_VCF_FILES/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"

    # # TRUTH BED
    # os.environ['TRUTH_BED_FILES_DIR']    = f"{os.environ['INPUT_DIR']}/TRUTH_BED_FILES"
    # os.environ['TRUTH_BED_FILE']         = f"{os.environ['INPUT_DIR']}/TRUTH_BED_FILES/HG002_GRCh38_1_22_v4.2.1_benchmark.bed"

    # # REFERENCE 
    # os.environ['REF_FILES_DIR']          = f"{os.environ['INPUT_DIR']}/REF_FILES"
    # os.environ['REF_FILE']               = f"{os.environ['INPUT_DIR']}/REF_FILES/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz"

    # # TRAIN TEST IMAGES
    # os.environ['TRAIN_OUTPUT_DIR']       = f"{os.environ['INPUT_DIR']}/PEPPER_TRAIN_IMAGES"
    # os.environ['TEST_OUTPUT_DIR']        = f"{os.environ['INPUT_DIR']}/PEPPER_TEST_IMAGES"

    # # MODELS
    os.environ['DORADO_MODELS']          = f"{os.environ['BASE_DIR']}/dorado_models"
    # os.environ['MODEL_OUTPUT_DIR']       = f"{os.environ['OUTPUT_DIR']}/trained_models"
    # os.environ['MODEL']                  = f"{os.environ['MODEL_OUTPUT_DIR']}/trained_models_07112024_195240/PEPPER_VARIANT_STEP_20000_checkpoint.pkl"

    # # MODEL EVALUATION
    # os.environ['EVAL_OUTPUT_DIR']        = f"{os.environ['OUTPUT_DIR']}/pepper_output/pepper_eval_model"
    # os.environ['VCF_FILE']               = f"{os.environ['EVAL_OUTPUT_DIR']}/PEPPER_VARIANT_FULL.vcf.gz"
    # os.environ['HAPPY_OUTPUT_DIR']       = f"{os.environ['OUTPUT_DIR']}/happy_outputs"
    # os.environ['HAPPY_OUTPUT_FILE']      = f"{os.environ['HAPPY_OUTPUT_DIR']}/HG002_pepper_model_07112024_195240"

    # # REPLACE MODEL
    # os.environ['PEPPER_OUTPUT_DIR']      = f"{os.environ['OUTPUT_DIR']}/pepper_deepvariant_output"

    # # MARGIN
    # os.environ['MARGIN_OUTPUT_DIR']      = f"{os.environ['OUTPUT_DIR']}/margin_output"
    # os.environ['MARGIN_OUTPUT_PREFIX']   = f"{os.environ['MARGIN_OUTPUT_DIR']}/{sample_file_name}"
    # os.environ['MARGIN_PARAMETER_FILE']  = "allParams.phase_vcf.ont.json"

    # # DEEPVARIANT
    # os.environ['DEEPVARIANT_OUTPUT_DIR'] = f"{os.environ['OUTPUT_DIR']}/deepvariant_output"
    # os.environ['DEEPVARIANT_VCF_FILE']   = f"{os.environ['DEEPVARIANT_OUTPUT_DIR']}/deepvariant.vcf.gz"
    # os.environ['DEEPVARIANT_GVCF_FILE']  = f"{os.environ['DEEPVARIANT_OUTPUT_DIR']}/deepvariant.g.vcf.gz"

    # SNIFFLES
    os.environ['SNIFFLES_OUTPUT_DIR']    = f"{os.environ['OUTPUT_DIR']}/sniffles_output"
    os.environ['SNIFFLES_VCF_FILE']      = f"{os.environ['SNIFFLES_OUTPUT_DIR']}/sniffles.vcf.gz"

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

def convert_pod5_to_bam() -> None:
    """Convert POD5 file to BAM file using dorado basecaller.
    
    Parameters:
        -x : cpu only
    """
    command = (
        "time dorado basecaller "
        "-x cpu "
        f"{os.environ['DORADO_MODELS']}/dna_r10.4.1_e8.2_400bps_fast@v4.2.0 "
        f"{os.environ['POD5_FILE']} > " 
        f"{os.environ['BAM_FILE']}"
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
        "time dorado basecaller "
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
        "time dorado basecaller "
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
        "time seqtk seq "
        f"-a {os.environ['FASTQ_FILE']} > " 
        f"{os.environ['FASTA_FILE']}"
    )
    run_command(command)

def create_fasta_index_file() -> None:
    """Create FASTA index file."""
    command = (
        "time samtools faidx "
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
        "time samtools faidx "
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
        "time dorado aligner "
        f"-t {os.environ['THREADS']} "
        f"   {os.environ['REF_FILE']} "
        f"   {os.environ['FASTA_FILE']} > "
        f"   {os.environ['BAM_FILE']}" 
    )
    run_command(command)

def create_bam_index_file() -> None:
    """Create bam index file."""
    command = (
        "time samtools index "
        f"-@ {os.environ['THREADS']} "
        f"   {os.environ['BAM_FILE']}"
    )
    run_command(command)

def sort_bam_file() -> None:
    """Sort bam file."""
    command = (
        "time samtools sort "
        f"--threads {os.environ['THREADS']} "
        f"          {os.environ['BAM_FILE']} "
        f"-o        {os.environ['BAM_SORTED_FILE']}"
    )
    run_command(command)

def create_sorted_bam_index_file() -> None:
    """Create sorted bam index file."""
    command = (
        "time samtools index "
        f"-@ {os.environ['THREADS']} "
        f"   {os.environ['BAM_SORTED_FILE']}"
    )
    run_command(command)

def create_truth_vcf_index_file() -> None:
    """Create truth vcf index file."""

    # delete existing vcf index files with suffix .tbi
    truth_vcf_files_dir = os.environ.get('TRUTH_VCF_FILES_DIR')

    if truth_vcf_files_dir:
        for filename in os.listdir(truth_vcf_files_dir):
            if filename.endswith('.tbi'):
                file_path = os.path.join(truth_vcf_files_dir, filename)
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)

    command = (
        "time tabix "
        f"{os.environ['TRUTH_VCF_FILE']}"
    )
    run_command(command)

def make_train_images(region: str, downsample_ratio: float, precision: float) -> None:
    """
    Run make_train_images command with specified parameters.

    Parameters:
        -b  : Input sorted BAM file.
        -f  : Input reference fasta file.
        -tv : Input truth VCF file.
        -r  : Regions (chr1-19 for training).
        -rb : GIAB high-confidence regions bed file.
        -o  : Output directory.
        -d  : Downsample fraction.
        -p  : Probability of a homozygous site being selected for training. This is to maintain class balance.
        --ont_r9_guppy5_sup : Use ONT presets for candidate finding.
    """

    # delete existing training image files
    train_image_files_dir = os.environ.get('TRAIN_OUTPUT_DIR')

    if train_image_files_dir:
        for filename in os.listdir(train_image_files_dir):
            file_path = os.path.join(train_image_files_dir, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)

    command = (
        "time pepper_variant_train make_train_images "
        f"-b   {os.environ['BAM_SORTED_FILE']} "
        f"-f   {os.environ['REF_FILE']} "
        f"-tv  {os.environ['TRUTH_VCF_FILE']} "
        f"-rb  {os.environ['TRUTH_BED_FILE']} "
        f"-o   {os.environ['TRAIN_OUTPUT_DIR']} "
        f"-t   {os.environ['THREADS']} "
        f"-r   {region} "
        f"-d   {downsample_ratio} "
        f"-p   {precision} "
        "--ont_r9_guppy5_sup"
    )
    run_command(command)

def make_test_images(region: str, downsample_ratio: float, precision: float) -> None:
    """
    Run make_train_images command with specified parameters.

    Parameters:
        -b  : Input sorted BAM file.
        -f  : Input reference fasta file.
        -tv : Input truth VCF file.
        -r  : Regions (chr20 for test).
        -rb : GIAB high-confidence regions bed file.
        -o  : Output directory.
        -d  : Downsample fraction.
        -p  : Probability of a homozygous site being selected for training. This is to maintain class balance.
        --ont_r9_guppy5_sup : Use ONT presets for candidate finding.
    """

    # delete existing testing image files
    test_image_files_dir = os.environ.get('TEST_OUTPUT_DIR')

    if test_image_files_dir:
        for filename in os.listdir(test_image_files_dir):
            file_path = os.path.join(test_image_files_dir, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)

    command = (
        "time pepper_variant_train make_train_images "
        f"-b   {os.environ['BAM_SORTED_FILE']} "
        f"-f   {os.environ['REF_FILE']} "
        f"-tv  {os.environ['TRUTH_VCF_FILE']} "
        f"-rb  {os.environ['TRUTH_BED_FILE']} "
        f"-o   {os.environ['TEST_OUTPUT_DIR']} "
        f"-t   {os.environ['THREADS']} "
        f"-r   {region} "
        f"-d   {downsample_ratio} "
        f"-p   {precision} "
        "--ont_r9_guppy5_sup"
    )
    run_command(command)

def train_model() -> None:
    """
    Run train_model command with specified parameters. Requires Nvidia GPU's.

    Must transfer training and test images from ubuntu development machine to AWS S3 bucket and then to AWS EC2.

    Parameters:
        -train : Path to train images directory.
        -test  : Path to test images directory.
        -o     : Output directory where models will be saved.
        -bs    : Batch size during training.
        -lr    : Learning rate value.
        -wd    : Weight decay value.
        -s     : Step size.
        -e     : Total epochs.
        -w     : Number of workers.
        --test_batch_size : Batch size during test (does not affect training).
        -g     : Set to use Nvidia GPUs.

    AWS:
        g4dn.xlarge, vCPU = 4, RAM = 16 MiB, Disk = SSD, Net bandwidth up to 25 Gbit / sec.

    Notes:
        This step may need high-speed bandwith to RAM.  And may need lots of RAM.
    """
    command = (
        "time pepper_variant_train train_model "
        f"-train {os.environ['TRAIN_OUTPUT_DIR']} "
        f"-test  {os.environ['TEST_OUTPUT_DIR']} "
        f"-o     {os.environ['MODEL_OUTPUT_DIR']} "
        "-bs     128 "
        "-lr     0.001 "
        "-wd     0.00001 "
        "-s      10000 "
        "-e      10 " # default 1000
        "-w      0 "
        "--test_batch_size 512 "
        "-g"
    )
    run_command(command)

def call_variant() -> None:
    """
    Run the call_variant command with the specified parameters.

    Parameters:
        -b  : Input BAM file.
        -f  : Input reference fasta file.
        -m  : Input model file.
        -o  : Output directory.
        -s  : Sample name.
        -t  : Number of threads.
        -w  : Window size.
        -bs : Batch size.
        -r  : Regions (chr20 for test).
        --ont_r9_guppy5_sup : Use ONT presets for candidate finding.
    """
    command = (
        "time pepper_variant call_variant "
        f"-b  {os.environ['BAM_FILE']} "
        f"-f  {os.environ['REF_FILE']} "
        f"-m  {os.environ['MODEL']} "
        f"-o  {os.environ['EVAL_OUTPUT_DIR']} "
        f"-t  {os.environ['THREADS']} "
        "-s   HG003 "
        "-w   0 "
        "-bs  512 "
        "-r   chr20 "
        "--ont_r9_guppy5_sup"
    )
    run_command(command)

def run_benchmark() -> None:
    """
    Run hap.py evaluation with specified parameters.
    
    Parameters:
        -f          : Truth BED file.
        -r          : Reference fasta file.
        -o          : Output directory.
        -l          : Regions (chr20 for test).
        --pass-only : Only evaluate passing variants.
        --no-roc    : Do not generate ROC curves.
        --no-json   : Do not generate JSON output.
        --engine    : Evaluation engine.
        --threads   : Number of threads.

    Ignore this warning message.  The reference file is set explicitely with the "-r" parameter.
    WARNING  No reference file found at default locations. You can set the environment variable 'HGREF' or 'HG19' to point to a suitable Fasta file.
    """
    command = (
        "time docker run "
        "-v    /home/seqcenter/PepperPipeline:/home/seqcenter/PepperPipeline "
        "-v    /data:/data "
        "-v    /data2:/data2 "
        "jmcdani20/hap.py:v0.3.12 "
        "/opt/hap.py/bin/hap.py "
        f"-f   {os.environ['TRUTH_BED_FILE']} "
        f"-r   {os.environ['REF_FILE']} "
        f"     {os.environ['TRUTH_VCF_FILE']} "
        f"     {os.environ['VCF_FILE']} "
        f"-o   {os.environ['HAPPY_OUTPUT_FILE']} "
        "-l    chr20 "
        "--pass-only "
        "--no-roc "
        "--no-json "
        "--engine=vcfeval "
        f"--threads={os.environ['THREADS']}"
    )
    run_command(command)
    
def run_replace_model() -> None:
    """
    Run the run_pepper_margin_deepvariant command with the specified parameters.
    
    Parameters: 
        -b  : Input BAM file.
        -f  : Input reference fasta file.
        -o  : Output directory.docs/html/
        -t  : Number of threads.
    """
    command = (
        "time docker run "
        f"-v   {os.environ['BASE_DIR']}:{os.environ['BASE_DIR']} "
        f"-v   {os.environ['INPUT_DIR']}:{os.environ['INPUT_DIR']} "
        f"-v   {os.environ['OUTPUT_DIR']}:{os.environ['OUTPUT_DIR']} "
        "kishwars/pepper_deepvariant:r0.8 "
        "run_pepper_margin_deepvariant call_variant "
        f"-b   {os.environ['BAM_FILE']} "
        f"-f   {os.environ['REF_FILE']} "
        f"-o   {os.environ['PEPPER_OUTPUT_DIR']} "
        f"-t   {os.environ['THREADS']} "
        f"--pepper_model {os.environ['MODEL']} "
        "--ont_r9_guppy5_sup"
    )
    run_command(command)

def run_margin() -> None:
    """
    Perform haplotyping with margin.
    
    Parameters:
        -o : Output prefix name.
        -t : Number of threads.
    """

    # we have to modify the paths to map the ubuntu host paths to the docker paths.  Should fix this eventually by replacing the docker file.
    base_path      = os.environ['BASE_DIR']
    bam_file       = os.environ['BAM_FILE'].replace(base_path, "")
    ref_file       = os.environ['REF_FILE'].replace(base_path, "")
    vcf_file       = os.environ['VCF_FILE'].replace(base_path, "")
    parameter_file = os.environ['MARGIN_PARAMETER_FILE']
    output_dir     = os.environ['MARGIN_OUTPUT_DIR'].replace(base_path, "")

    command = (
        "time docker run "
        "-v /home/seqcenter/PepperPipeline/data:/data "
        "pepper_deepvariant "
        "margin phase "
        f"{bam_file} "
        f"{ref_file} "
        # f"{vcf_file} "
        f"HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz "
        f"/opt/margin_dir/params/phase/{parameter_file} "
        f"-o {output_dir} "
        f"-t {os.environ['THREADS']}"
    )
    run_command(command)

def run_deepvariant() -> None:
    """
    Perform variant calling.
    
    Parameters:
        --model_type       : Replace this string with exactly one of the following [WGS,WES,PACBIO,ONT_R104,HYBRID_PACBIO_ILLUMINA].
        --ref              : Reference genome.
        --reads            : Input BAM file.
        --output_vcf       : Output VCF file. 
        --output_gvcf      : Output gVCF file
        --num_shards       : Number of threads.
        --logging_dir      : Optional. This saves the log output for each stage separately.
        --haploid_contigs  : Optional. Heterozygous variants in these contigs will be re-genotyped as the most likely of reference or homozygous alternates. For a sample with karyotype XY, it should be set to "chrX,chrY" for GRCh38 and "X,Y" for GRCh37. For a sample with karyotype XX, this should not be used.
        --par_regions_bed= : Optional. If --haploid_contigs is set, then this can be used to provide PAR regions to be excluded from genotype adjustment. Download links to this files are available in this page.
        --par_regions_bed="/input/GRCh3X_par.bed" " : Optional. If --haploid_contigs is set, then this can be used to provide PAR regions to be excluded from genotype adjustment. Download links to this files are available in this page.
        --dry_run          : Default is false. If set to true, commands will be printed out but not executed.
    """

    # we have to modify the paths to map the ubuntu host paths to the docker paths.  Should fix this eventually by replacing the docker file.
    train_input_path  = os.environ['INPUT_DIR']
    train_output_path = os.environ['OUTPUT_DIR']
    ref_file          = os.environ['REF_FILE'].replace(train_input_path, "")
    bam_sorted_file   = os.environ['BAM_SORTED_FILE'].replace(train_input_path, "")
    vcf_file          = os.environ['DEEPVARIANT_VCF_FILE'].replace(train_output_path, "")
    gvcf_file         = os.environ['DEEPVARIANT_GVCF_FILE'].replace(train_output_path, "")

    # BIN_VERSION="1.6.0"
    command = (
        "time docker run "
        f"-v {os.environ['INPUT_DIR']}:/input "             
        f"-v {os.environ['OUTPUT_DIR']}:/output "   
        "google/deepvariant:1.6.0 "                       # BIN_VERSION
        "/opt/deepvariant/bin/run_deepvariant "
        "--model_type=WGS "                               # Replace this string with exactly one of the following [WGS,WES,PACBIO,ONT_R104,HYBRID_PACBIO_ILLUMINA]
        # f"--ref=/input/{os.environ['REF_FILE']} "
        # f"--reads=/input/{os.environ['BAM_FILE']} "
        # f"--output_vcf=/output/{os.environ['DEEPVARIANT_VCF_FILE']} "
        # f"--output_gvcf=/output/{os.environ['DEEPVARIANT_GVCF_FILE']} "
        f"--ref=/input{ref_file} "
        f"--reads=/input{bam_sorted_file} "
        f"--output_vcf=/output{vcf_file} "
        f"--output_gvcf=/output{gvcf_file} "              # Optional. genomic VCF file.
        f"--num_shards={os.environ['THREADS']} "
        # "--logging_dir     = /output/logs "             # Optional. This saves the log output for each stage separately.
        # "--haploid_contigs = "chrX,chrY" "              # Optional. Heterozygous variants in these contigs will be re-genotyped as the most likely of reference or homozygous alternates. For a sample with karyotype XY, it should be set to "chrX,chrY" for GRCh38 and "X,Y" for GRCh37. For a sample with karyotype XX, this should not be used.
        # "--par_regions_bed = "/input/GRCh3X_par.bed" "  # Optional. If --haploid_contigs is set, then this can be used to provide PAR regions to be excluded from genotype adjustment. Download links to this files are available in this page.
        "--dry_run = false"                               # Default is false. If set to true, commands will be printed out but not executed.
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
        "time sniffles "
        f"-i {os.environ['BAM_SORTED_FILE']} "
        f"-v {os.environ['SNIFFLES_VCF_FILE']} "
        f"-t {os.environ['THREADS']} "
        "--allow-overwrite"
    )
    run_command(command)

# main entry point
if __name__ == "__main__":
    try:
        logger.info("Start PepperPipeline...")

        # set env vars
        logger.info("Set environment variables...")
        set_environment_variables()

        # # delete all files in subdirectories. Use with caution.
        # # runtime: 1 min.
        # logger.info("Clear files...")
        # clear_files()

        # # Add new mappings
        # add_mapping("test",       "data/input/test/input3.txt",       "data/output/test/output3.txt")
        # add_mapping("production", "data/input/production/input3.txt", "data/output/production/output3.txt")
        
        # # Get output file for a given input file
        # output_file = get_output_file("test", "data/input/test/{sample_file_name}.pod5")
        # print(f"Output file for 'data/output/test/{sample_file_name}.pod5': {output_file}")
        # output_file = get_output_file("production", "data/input/production/input3.txt")
        # print(f"Output file for 'data/input/production/input3.txt': {output_file}")

        # # TODO: replace with real datasets
        # # download alignment files
        # logger.info("Download alignment files...")
        # try:
        #     wget.download("https://storage.googleapis.com/pepper-deepvariant-public/pepper_deepvariant_training/HG002_guppy422_2_GRCh38_no_alt.bam",     out=os.environ['BAM_FILES_DIR'])
        #     wget.download("https://storage.googleapis.com/pepper-deepvariant-public/pepper_deepvariant_training/HG002_guppy422_2_GRCh38_no_alt.bam.bai", out=os.environ['BAM_FILES_DIR'])
        # except Exception as e:
        #     logger.error(f"ERROR: an error occurred while downloading alignment files: {e}")
        #     sys.exit(1)

        # # TODO: replace with real datasets
        # # download reference files
        # logger.info("Download reference files...")
        # try:
        #     wget.download("https://storage.googleapis.com/kishwar-helen/variant_calling_data/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",     out=os.environ['REF_FILES_DIR'])
        #     wget.download("https://storage.googleapis.com/kishwar-helen/variant_calling_data/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai", out=os.environ['REF_FILES_DIR'])
        # except Exception as e:
        #     logger.error(f"ERROR: an error occurred while downloading reference files: {e}")
        #     sys.exit(1)

        # # TODO: replace with real datasets
        # # download GIAB truth files
        # logger.info("Download GIAB truth files...")
        # try:
        #     wget.download("https://storage.googleapis.com/pepper-deepvariant-public/pepper_deepvariant_training/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",     out=os.environ['TRUTH_VCF_FILES_DIR'])
        #     wget.download("https://storage.googleapis.com/pepper-deepvariant-public/pepper_deepvariant_training/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi", out=os.environ['TRUTH_VCF_FILES_DIR'])
        #     wget.download("https://storage.googleapis.com/pepper-deepvariant-public/pepper_deepvariant_training/HG002_GRCh38_1_22_v4.2.1_benchmark.bed",        out=os.environ['TRUTH_BED_FILES_DIR'])
        # except Exception as e:
        #     logger.error(f"ERROR: an error occurred while downloading GIAB truth files: {e}")
        #     sys.exit(1)

        # # 1. convert pod5 file to bam file
        # # runtime: 5 min.
        # logger.info("Convert pod5 to bam...")
        # convert_pod5_to_bam()

        # # convert fast5 file to bam file
        # # runtime: 5 min.
        # logger.info("Convert fast5 to bam...")
        # convert_fast5_to_bam()

        # # 1. convert pod5 file to fastq file
        # # runtime: 5 min.
        # logger.info("Convert pod5 to fastq...")
        # convert_pod5_to_fastq()

        # # 2. convert fastq file to fasta file
        # # runtime: 5 min.
        # logger.info("Convert fastq to fasta...")
        # convert_fastq_to_fasta()

        # # 3. create fasta index file
        # # runtime: 3 min.
        # logger.info("Create fasta index file...")
        # create_fasta_index_file()

        # # 4. create reference genome index file
        # # runtime: 3 min.(peppervenv) seqcenter@seqcenter-MS-7C75:~/PepperPipelin
        # logger.info("Create reference genome index file...")
        # create_ref_genome_index_file()

        # # 5. align fasta file to reference genome
        # # runtime: 3 min.
        # logger.info("Align fasta to reference genome...")
        # align_fasta_to_reference()

        # # 6. create bam index file
        # # runtime: 1 min.
        # logger.info("Create bam index file...")
        # create_bam_index_file()

        # 2. sort bam file
        # runtime: 1 min.
        logger.info("Sort bam file...")
        sort_bam_file()

        # # 8. create sorted bam index file
        # # runtime: 1 min.
        # logger.info("Create sorted bam index file...")
        # create_sorted_bam_index_file()

        # # 9. create truth vcf index file
        # # runtime: 1 min.
        # logger.info("Create truth vcf index file...")
        # create_truth_vcf_index_file()

        # # 10. generate training images
        # # runtime: 158 min. (2.63 hrs.)
        # logger.info("Generate training images...")
        # make_train_images(
        #     region           = 'chr1-19',
        #     downsample_ratio = 1.0,
        #     precision        = 0.3
        # )

        # # 11. generate testing images
        # # runtime: 3 min.
        # logger.info("Generate testing images...")
        # make_test_images(
        #     region           = 'chr20',
        #     downsample_ratio = 1.0,
        #     precision        = 1.0
        # )

        # # 12. train model with training images. requires nvidia gpu's    #     f"   {os.environ['MARGIN_PARAMETER_FILE']} "
        # # runtime: 
        # logger.info("Train model with training images...")
        # train_model()

        # # 13. generate variant calls
        # # runtime: 20 min.
        # logger.info("Generate variant calls...")
        # call_variant()

        # # 14. benchmark the model
        # # runtime: 3 min.
        # logger.info("Benchmark the model...")
        # run_benchmark()

        # # 15. perform haplotyping
        # # runtime: 3 min.
        # logger.info("Perform haplotyping...")
        # run_margin()

        # # 16. perform variant calling 
        # # runtime: 3 min.
        # logger.info("Perform variant calling...")
        # run_deepvariant()

        # 3. perform structural variant calling
        # runtime: 10 min.
        logger.info("Perform structural variant calling...")
        run_sniffles()

        # # replace default model with trained model
        # # runtime: 860 min. (14 hrs.)
        # logger.info("Replace default model with trained model...")
        # run_replace_model()

        logger.info("Finish PepperPipeline...")

    except Exception as e:
        logger.error(f"ERROR: an error occurred in main: {e}")
        sys.exit(1)
