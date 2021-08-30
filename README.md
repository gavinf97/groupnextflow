# Group Nextflow Assignment

## Assignment Description
Produce a working Hi-C nextflow script based off of the NF-core Hi-C pipeline steps and also include a singularity container for reproducibility.
https://github.com/nf-core/hic

****Note: all code is strongly based off of the nf-core Hi-C script linked****

## Description of files for runnning Hi-C

### 1. Nextflow Script
'hi-c.nf'
Runs Hi-C analsyis in conjunction with Singularity container. Sample input FASTQ files present for test run.

### Nextflow Dependency Files
#### 1.1 FASTQ Input Samples
For use in the Hi-C pipeline to produce an output
Test dataset reference source: https://github.com/nf-core/test-datasets/tree/hic

#### 1.2 Alignment scripts
align.nf
-Run “align.nf” to do Step 1 of the pipeline (mapping using a two steps strategy to rescue reads spanning the ligation sites). 

#### 1.3 Python Merging Script
mergeSAM.py (any other required python scripts) 
-"mergeSAM.py" merges the SAM files (which have undergone the 2-step alignment) into one paired-end BAM file – the final output of Step 1. 
This script and others located here: https://github.com/nf-core/hic/blob/master/bin.

#### 1.4 NF Step one Output File Dependency
Final output file of step 1 is ${sample}_bwt2pairs.bam (SRR4292758_00_bwt2pairs.bam).

### 2. Singularity Container Directories and Descriptions
Base availability: Dockerfile and environment YAML: https://github.com/nf-core/hic


****Modified:****
Index files placed in '/results/index' directory, all other files are placed in '/results/align' directory.
'/data' directory: contains the raw fastq.gz files (SRR4292758_00_R1.fastq.gz, SRR4292758_00_R2.fastq.gz)
'/reference' directory: contains the reference genome (W303_SGD_2015_JRIU00000000.fsa). N.B. - reference used to be .fsa.txt


