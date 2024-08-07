# Sniffles2
A fast structural variant caller for long-read sequencing, Sniffles2 accurately detect SVs on germline, somatic and population-level for PacBio and Oxford Nanopore read data.

## Quick Start: Germline SV calling using Sniffles2
To call SVs from long read alignments (PacBio / ONT), you can use:

`sniffles -i mapped_input.bam -v output.vcf`

For improved calling in repetitive regions, Sniffles2 accepts a tandem repeat annotations file using the option `--tandem-repeats annotations.bed`. Sniffles2 compatible tandem repeat annotations for human references can be downloaded from the [annotations/ folder](https://github.com/fritzsedlazeck/Sniffles/tree/master/annotations).

(see sniffles --help or below for full usage information).

## Installation
Install Sniffles2 using pip:
```
pip install --no-cache-dir sniffles`
```


## Requirements
* Python >= 3.10
* pysam >= 0.21.0
* edlib >=1.3.9
* psutil>=5.9.4

Notes:   In production mode, this application requires Nvidia GPU's (H100, A100, V100).

#### Tested on:
* python==3.10.12
* pysam==0.21.0

<br>

## AWS (Amazon Web Services)  
EC2  
PepperPipeline is hosted on AWS EC2 instances that support Nvidia GPU's.  Here are some commonly used instance types.  In all cases the architecture is x86_64.  

| Instance type | # vCPU's | Clock speed (GHz) | CPU Memory (GiB)  | Storage (GB) | Storage type | Network Performance (Gbit/sec.) | GPU name | # GPU's | GPU memory (GiB) | Price (USD/hr.) |  
| :-----------: | -------: | ----------------: | ----------------: | -----------: | -----------: | ------------------------------: | -------: | ------: | ---------------: | --------------: |  
| g4dn.xlarge   |     4    |      2.5          |           16      |    125       | SSD          | Up to 25                        | T4       | 1       | 16               |  0.71           |  
|   g5.xlarge   |     4    |      3.3          |           16      |    250       | SSD          | Up to 10                        | A10G     | 1       | 24               |  1.01           |  
|   g6.xlarge   |     4    |      3.4          |           16      |    250       | SSD          | Up to 10                        | L4       | 1       | 22               |  0.80           |  
|   p3.2xlarge  |     8    |      2.7          |           61      |    ---       | ---          | Up to 10                        | V100     | 1       | 16               |  3.06           |  
|  p4d.24xlarge |    96    |      3.0          |        1,152      |  8,000       | SSD          | 4x 100                          | A100     | 8       | 40               | 32.77           |  
|   p5.48xlarge |   192    |      3.6          |        2,048      | 30,400       | SSD          | 3,200                           | H100     | 8       | 80               | 98.32           |  

<br>

## Citation
Please cite our paper at:
Sniffles v2: 
https://www.nature.com/articles/s41587-023-02024-y

and 
Sniffles v1:
https://www.nature.com/articles/s41592-018-0001-7

## Use-Cases / Modes

### A. General (all Modes)
* To output deletion (DEL SV) sequences, the reference genome (.fasta) must be specified using e.g. `--reference reference.fasta`.
* Sniffles2 supports optionally specifying tandem repeat region annotations (.bed), which can improve calling in these regions `--tandem-repeats annotations.bed`. Sniffles2 compatible tandem repeat annotations for human references can be found in the [annotations/ folder](https://github.com/fritzsedlazeck/Sniffles/tree/master/annotations).
* Sniffles2 is fully parallelized and uses 4 threads by default. This value can be adapted using e.g. `--threads 4` as option. Memory requirements will increase with the number of threads used.
* To output read names in SNF and VCF files, the `--output-rnames` option is required.

### B. Multi-Sample SV Calling (Trios, Populations)
Multi-sample SV calling using Sniffles2 population mode works in two steps:

1. Call SV candidates and create an associated .snf file for each sample: `sniffles --input sample1.bam --snf sample1.snf`
2. Combined calling using multiple .snf files into a single .vcf: `sniffles --input sample1.snf sample2.snf ... sampleN.snf --vcf multisample.vcf`

Alternatively, for step 2. you can supply a .tsv file, containing a list of .snf files, and custom sample ids in an optional second column (one sample per line), .e.g.:
2. Combined calling using a .tsv as sample list: `sniffles --input snf_files_list.tsv --vcf multisample.vcf`

### C. Mosaic SV Calling (Non-germline or somatic SVs)
To call mosaic SVs, the *--mosaic* option should be added, i.e.:

`sniffles --input mapped_input.bam --vcf output.vcf --mosaic`

### D. Genotyping a known set of SVs (Force Calling)
Example command, to determine the genotype of each SV in *input_known_svs.vcf* for *sample.bam* and write the re-genotyped SVs to *output_genotypes.vcf*:

`sniffles --input sample.bam --genotype-vcf input_known_svs.vcf --vcf output_genotypes.vcf`

## Quick Tips

### Input / Output
* .bam or .cram files containing long read alignments (i.e. from minimap2 or ngmlr) are supported as input
* .vcf.gz (bgzipped+tabix indexed) output is supported
* Simultaneous output of both .vcf and .snf file (for multi-sample calling) is supported

## Companion apps
* We have developed a plotting tools for Sniffles2: [https://github.com/farhangus/sniffle2_plot](https://github.com/farhangus/sniffle2_plot)
* We also provide VCF and scripts used for the manuscript [https://github.com/smolkmo/Sniffles2-Supplement](https://github.com/smolkmo/Sniffles2-Supplement) 

## Supplementary tables
[https://github.com/smolkmo/Sniffles2-Supplement/blob/main/Supplemetary%20tables.xlsx](https://github.com/smolkmo/Sniffles2-Supplement/blob/main/Supplemetary%20tables.xlsx)
