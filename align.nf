#!/bin/bash/env nextflow

reads = "$baseDir/data/*_R{1,2}.fastq.gz"
reference = "$baseDir/reference/*.fsa.txt"
outDir = "$baseDir/results"
project_dir = projectDir

Channel
	.fromPath(reads)
	.map{file -> [file.simpleName, file]}
	.set{reads_ch}

Channel
	.fromPath(reference)
	.map{file -> [file.simpleName, file]}
	.set{reference_ch}


process create_bowtie2_index {
	publishDir path: "$outDir/index", mode: "copy"

	input:
	tuple val(base), file(ref) from reference_ch

	output:
	tuple val(base), file("*.bt2") into index_ch
	tuple val(base), file("*.bt2") into index_trimmed_ch

	script:
	"""
	bowtie2-build ${ref} ${base}
	"""
}


process bowtie2_end_to_end {
	publishDir path: "$outDir/align", mode: "copy"

	input:
	each(reads) from reads_ch
	tuple val(base), file(index) from index_ch

	output:
	tuple val(prefix), file("${reads[0]}.sam") into end_to_end_sam_ch
	tuple val(prefix), file(name) into unmapped_ch

	script:
	prefix = reads[0]
	name = "${prefix}_unmap.fastq"

        """
        bowtie2 -x ${base} \\
        --un ${reads[0]}_unmap.fastq \\
        -U ${reads[1]} \\
        -S ${reads[0]}.sam \\
        --very-sensitive \\
        -L 30 \\
        --score-min L,-0.6,-0.2 \\
        --end-to-end \\
        --reorder
        """
}


process trim_reads {
	publishDir path: "$outDir/align", mode: "copy"

        input:
        each(reads) from unmapped_ch

        output:
        tuple val(file), file(name) into trimmed_ch

        script:
	file = reads[0]
	name = "${file}_unmap.fastq.trimmed"

        """
        homerTools trim -3 AAGCTAGCTT \\
        	-mis 0 \\
        	-matchStart 20 \\
        	-min 20 \\
        	${reads[1]}

	mv ${reads[1]}.trimmed ${name}
        """
}


process bowtie2_on_trimmed_reads {
	publishDir path: "$outDir/align", mode: "copy"

	input:
	each(reads) from trimmed_ch
	tuple val(name), file(index) from index_trimmed_ch

        output:
        tuple val(prefix), file("${reads[0]}.trimmed.sam") into trimmed_sam_ch

        script:
	prefix = reads[0]

        """
        bowtie2 -x ${name} \\
        -U ${reads[1]} \\
        -S ${reads[0]}.trimmed.sam \\
        --very-sensitive \\
        -L 20 \\
        --score-min L,-0.6,-0.2 \\
        --end-to-end \\
        --reorder
        """
}


process merge_mapping_steps{
	publishDir path: "$outDir/align", mode: "copy"

	input:
	tuple val(name), file(bam1), file(bam2) from end_to_end_sam_ch.join(trimmed_sam_ch)

	output:
	tuple val(prefix), file("${name}_bwt2merged.sam") into bwt2_merged_sam

	script:
	prefix = name.toString() - ~/(_R1|_R2|_val_1|_val_2|_1|_2)$/

	"""
	samtools merge -f ${name}_bwt2merged.sam \\
	${bam1} ${bam2}

	samtools sort -n -T /tmp/ \\
	-o ${name}_bwt2merged.sorted.sam \\
	${name}_bwt2merged.sam

	mv ${name}_bwt2merged.sorted.sam ${name}_bwt2merged.sam
	"""
}


process combine_mapped_files{
	publishDir path: "$outDir/align", mode: "copy"

	input:
	tuple val(sample), file(aligned_sam) from bwt2_merged_sam.groupTuple()

	output:
	tuple val(sample), file("${sample}_bwt2pairs.bam") into paired_bam

	script:
	"""
	python $project_dir/mergeSAM.py \\
	-f ${aligned_sam[0]} \\
	-r ${aligned_sam[1]} \\
	-o ${sample}_bwt2pairs.bam
	"""
}
