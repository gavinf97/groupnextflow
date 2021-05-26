#!/bin/bash/env nextflow

/*
 * Written by Sarah
 */

reads = "$baseDir/data/*_R{1,2}.fastq.gz"
reference = "$baseDir/reference/*.fsa"
outDir = "$baseDir/results"
project_dir = projectDir

Channel
	.fromPath(reads)
	.map{file -> [file.simpleName, file]}
	.into{ reads_ch; reads_for_fastqc }

Channel
	.fromPath(reference)
	.map{file -> [file.simpleName, file]}
	.into{ reference_ch; fasta_for_resfrag_ch; fasta_for_chromsize }


process create_bowtie2_index {
	publishDir path: "$outDir/index", mode: "copy"

	input:
	tuple val(base), file(ref) from reference_ch

	output:
	tuple val(base), file("*.bt2") into index_ch

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
	tuple val(prefix), file("${name}.bam") into end_to_end_bam

	script:
	name = reads[0]
	prefix = name.toString() - ~/(_R1|_R2|_val_1|_val_2|_1|_2)$/

        """
        bowtie2 -x ${base} \\
        -U ${reads[1]} \\
        -S ${reads[0]}.bam \\
        --very-sensitive \\
        -L 30 \\
        --score-min L,-0.6,-0.2 \\
        --end-to-end \\
        --reorder
        """
}

process combine_mapped_files{
	publishDir path: "$outDir/align", mode: "copy",
		saveAs: {filename ->
		filename.indexOf(".pairstat") > 0 ? "$outDir/stats/$filename" : "$filename"}

	input:
	tuple val(sample), file(aligned_bam) from end_to_end_bam.groupTuple()

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

/*
 * Written by Gavin
 */

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


/*
 * Written by Aodán
 */

process make_chrom_size {
        publishDir path: "$outDir/chrom_size", mode: "copy"

        input:
        tuple val(base), file(fasta) from fasta_for_chromsize

        output:
        file("*.size") into chromosome_size, chromosome_size_cool

        script:
        """
	samtools faidx ${fasta}
	cut -f1,2 ${fasta}.fai > chrom.size
	"""
}



// Resolutions for contact maps
bin_size = '1000000,500000'
map_res = bin_size.tokenize(',')

process build_contact_maps {
	publishDir path: "$outDir/matrix/raw", mode: "copy"

	input:
	tuple val(sample), file(vpairs), val(mres) from all_valid_pairs.combine(map_res)
	file chrsize from chromosome_size.collect()

	output:
	file("*.matrix") into raw_maps
	file "*.bed"

	script:
	"""
	${project_dir}/build_matrix --matrix-format upper  \\
	--binsize ${mres} \\
	--chrsizes ${chrsize} \\
	--ifile ${vpairs} \\
	--oprefix ${sample}_${mres}
	"""
}


process run_ice{
	publishDir "$outDir/matrix/iced", mode: "copy"

	input:
	file(rmaps) from raw_maps

	output:
	file("*iced.matrix") into iced_maps

	script:
	prefix = rmaps.toString() - ~/(\.matrix)?$/
	"""
	ice --filter_low_counts_perc 0.02 \
	--results_filename ${prefix}_iced.matrix \
	--filter_high_counts_perc 0 \
	--max_iter 100 --eps 0.1 --remove-all-zeros-loci --output-bias 1 --verbose 1 ${rmaps}
	"""
}


/*
 * Written by Eric
 */

process merge_sample {
	publishDir "$outDir/mstats", mode: "copy"

	input:
	tuple val(prefix), file(fstat) from all_rsstat.groupTuple().concat(all_pairstat.groupTuple())

	output:
	file("*") into all_mstats

	script:
	sample = prefix.toString() - ~/(_R1|_R2|_val_1|_val_2|_1|_2)/
	if ( (fstat =~ /.pairstat/) ){ ext = "mpairstat" }
	if ( (fstat =~ /.RSstat/) ){ ext = "mRSstat" }

	"""
        mkdir -p $outDir/mstats/${sample}
        python $project_dir/merge_statfiles.py -f ${fstat} > ${prefix}.${ext}
	"""
}


process fastqc {
	input:
	tuple val(sample), file(reads) from reads_for_fastqc

	output:
	file "*_fastqc.{zip,html}" into fastqc_results

	script:
	"""
	fastqc -q ${reads}
	"""
}

process multiqc {
	publishDir "$outDir", mode: "copy"

	input:
	file ('*') from fastqc_results.collect().ifEmpty([])

	output:
	file "multiqc_report.html" into multiqc_report
	file "multiqc_data"

	script:
	"""
	multiqc .
	"""
}
