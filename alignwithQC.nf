#!/bin/bash/env nextflow

reads = "$baseDir/data/*_R{1,2}.fastq.gz"
reference = "$baseDir/reference/*.fsa"
outDir = "$baseDir/results"
project_dir = projectDir

Channel
	.fromPath(reads)
	.map{file -> [file.simpleName, file]}
	.set{reads_ch}

Channel
	.fromPath(reference)
	.map{file -> [file.simpleName, file]}
	.into{ reference_ch; fasta_for_resfrag_ch }


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


process get_restriction_fragments {
	publishDir path: "$outDir/index", mode: "copy"

	input:
	file(fasta) from fasta_for_resfrag_ch

	output:
	file("*.bed") into res_frag_ch

	script:
	"""
	python $project_dir/digest_genome.py \\
	-r A^AGCTT \\
	-o restriction_fragments.bed ${fasta[1]}
	"""
}


process bowtie2_end_to_end {
	publishDir path: "$outDir/align", mode: "copy"

	input:
	each(reads) from reads_ch
	tuple val(base), file(index) from index_ch

	output:
	tuple val(prefix), file("${reads[0]}.bam") into end_to_end_bam_ch
	tuple val(prefix), file(name) into unmapped_ch

	script:
	prefix = reads[0]
	name = "${prefix}_unmap.fastq"

        """
        bowtie2 -x ${base} \\
        --un ${reads[0]}_unmap.fastq \\
        -U ${reads[1]} \\
        -S ${reads[0]}.bam \\
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
	tuple val(prefix), file("${reads[0]}.trimmed.bam") into trimmed_bam_ch

	script:
	prefix = reads[0]

	"""
	bowtie2 -x ${name} \\
	-U ${reads[1]} \\
	-S ${reads[0]}.trimmed.bam \\
	--very-sensitive \\
	-L 20 \\
	--score-min L,-0.6,-0.2 \\
	--end-to-end \\
	--reorder
	"""
}


process merge_mapping_steps{
	publishDir path: "$outDir/align", mode: "copy",
		saveAs: {filename ->
		filename.indexOf(".mapstat") > 0 ? "$outDir/stats/$filename" : "$filename"}

	input:
	tuple val(prefix), file(bam1), file(bam2) from end_to_end_bam_ch.join(trimmed_bam_ch)

	output:
	tuple val(sample), file("${prefix}_bwt2merged.bam") into bwt2_merged_bam
	tuple val(oname), file("${prefix}.mapstat") into all_mapstat

	script:
	sample = prefix.toString() - ~/(_R1|_R2|_val_1|_val_2|_1|_2)$/
	oname = prefix.toString() - ~/(\.[0-9]+)$/
	tag = "$sample = $bam1 + $bam2"

	"""
	samtools merge -f ${prefix}_bwt2merged.bam \\
	${bam1} ${bam2}
	samtools sort -n -T /tmp/ \\
	-o ${prefix}_bwt2merged.sorted.bam \\
	${prefix}_bwt2merged.bam
	mv ${prefix}_bwt2merged.sorted.bam ${prefix}_bwt2merged.bam
	echo "## ${prefix}" >> ${prefix}.mapstat
	echo -n "total_${tag}\t" >> ${prefix}.mapstat
	samtools view -c ${prefix}_bwt2merged.bam >> ${prefix}.mapstat
	echo -n "mapped_${tag}\t" >> ${prefix}.mapstat
	samtools view -c -F 4 ${prefix}_bwt2merged.bam >> ${prefix}.mapstat
	echo -n "global_${tag}\t" >> ${prefix}.mapstat
	samtools view -c -F 4 ${bam1} >> ${prefix}.mapstat
	echo -n "local_${tag}\t"  >> ${prefix}.mapstat
	samtools view -c -F 4 ${bam2} >> ${prefix}.mapstat
	"""
}


process combine_mapped_files{
	publishDir path: "$outDir/align", mode: "copy",
		saveAs: {filename ->
		filename.indexOf(".pairstat") > 0 ? "$outDir/stats/$filename" : "$filename"}

	input:
	tuple val(sample), file(aligned_bam) from bwt2_merged_bam.groupTuple()

	output:
	tuple val(sample), file("${sample}_bwt2pairs.bam") into paired_bam
	tuple val(oname), file("*.pairstat") into all_pairstat

	script:
	oname = sample.toString() - ~/(\.[0-9]+)$/

	"""
	python $project_dir/mergeSAM.py \\
	-f ${aligned_bam[0]} \\
	-r ${aligned_bam[1]} \\
	-o ${sample}_bwt2pairs.bam \\
	--single --multi -t
	"""
}


 process get_valid_interaction{
	publishDir path: "$outDir/valid", mode: "copy",
		saveAs: {filename ->
		filename.indexOf(".RSstat") > 0 ? "$outDir/stats/$filename" : "$filename"}

	input:
	tuple val(sample), file(pe_bam) from paired_bam
	file(frag_file) from res_frag_ch.collect()

	output:
	tuple val(sample), file("*.validPairs") into valid_pairs
	tuple val(sample), file("*.validPairs") into valid_pairs_4cool
	tuple val(sample), file("*.RSstat") into all_rsstat

	script:
	"""
	python $project_dir/mapped_2hic_fragments.py \\
	-f ${frag_file} \\
	-r ${pe_bam} \\
	sort -T /tmp/ -k2,2V -k3,3n -k5,5V -k6,6n
	"""
}


process remove_duplicates {
	publishDir path: "$outDir/valid", mode: "copy",
		saveAs: {filename ->
		filename.indexOf(".mergestat") > 0 ? "$outDir/stats/$filename" : "$filename"}

	input:
	tuple val(sample), file(vpairs) from valid_pairs.groupTuple()

	output:
	tuple val(sample), file("*.allValidPairs") into all_valid_pairs
	tuple val(sample), file("*.allValidPairs") into all_valid_pairs_4cool
	file("*") into all_mergestat

	script:
	"""
	## Sort valid pairs and remove read pairs with same starts (i.e duplicated read pairs)
	sort -T /tmp/ -S 50% -k2,2V -k3,3n -k5,5V -k6,6n -m ${vpairs} | \
	awk -F"\\t" 'BEGIN{c1=0;c2=0;s1=0;s2=0}(c1!=\$2 || c2!=\$5 || s1!=\$3 || s2!=\$6){print;c1=\$2;c2=\$5;s1=\$3;s2=\$6}' >> ${sample}.allValidPairs
	echo -n "valid_interaction\t" >> ${sample}_allValidPairs.mergestat
	cat ${vpairs} | wc -l >> ${sample}_allValidPairs.mergestat
	echo -n "valid_interaction_rmdup\t" >> ${sample}_allValidPairs.mergestat
	cat ${sample}.allValidPairs | wc -l >> ${sample}_allValidPairs.mergestat
	## Count short range (<20000) vs long range contacts
	awk 'BEGIN{cis=0;trans=0;sr=0;lr=0} \$2 == \$5{cis=cis+1; d=\$6>\$3?\$6-\$3:\$3-\$6; if (d<=20000){sr=sr+1}else{lr=lr+1}} \$2!=\$5{trans=trans+1}END{print "trans_interaction\\t"trans"\\ncis_interaction\\t"cis"\\ncis_shortRange\\t"sr"\\ncis_longRange\\t"lr}' ${sample}.allValidPairs >> ${sample}_allValidPairs.mergestat
	"""
}

process fastqc {
    input:
    each(reads) from reads_ch

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}

process multiqc {
    input:
    file ('fastqc/*') from fastqc_results.collect().ifEmpty([])

    output:
    file "multiqc_report.html" into multiqc_report
    file "multiqc_data"

    script:
    """
    multiqc .
    """
}
