# groupnextflow
Hi C nextflow script and singularity container

Test dataset/reference source: https://github.com/nf-core/test-datasets/tree/hic

Within working directory: align.nf, mergeSAM.py, (any other required python scripts) and the subdirectories /data and /reference

mergeSAM.py merges the SAM files (which have undergone the 2-step alignment) into one paired-end BAM file – the final output of step 1. 
This script and others located here: https://github.com/nf-core/hic/blob/master/bin.

/data directory: contains the raw fastq.gz files (SRR4292758_00_R1.fastq.gz, SRR4292758_00_R2.fastq.gz)

/reference directory: contains the reference genome (W303_SGD_2015_JRIU00000000.fsa). N.B. - reference used to be .fsa.txt, I've changed the script too

Run “align.nf” to do step 1 of the pipeline (mapping using a two steps strategy to rescue reads spanning the ligation sites.) 
N.B. - I will be changing the trim_reads step so that other restriction enzymes can be used – this will be specified as a command line argument

Index files placed in /results/index, all other files are placed in /results/align. 
Final output file of step 1 is ${sample}_bwt2pairs.bam (SRR4292758_00_bwt2pairs.bam).
