#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.bam_directory_path = ''
params.bed_file_path = ''
params.fasta_file_path = ''
params.out_dir = ''
params.at_anchor = 25
params.gc_anchor = 75
params.failure_fold_change = 2
params.warning_fold_change = 1.5
params.failure_at = 1.5
params.failure_gc = 1.5
params.draw_trend = false
params.show_sample_names = true

// Check if the paths are not empty and the files exist
if (!params.bam_directory_path || !new File(params.bam_directory_path).exists()) {
    exit 1,  """
    Error: Missing or incorrect path for --bam_directory_path parameter.
    Please ensure the input directory for alignment BAM files exists.
    """
}
if (!params.bed_file_path || !new File(params.bed_file_path).exists()) {
    exit 1, """
    Error: Missing or incorrect path for --bed_file_path parameter.
    Please ensure the BED file for probes or genomic bins exists.
    """
}
if (!params.fasta_file_path || !new File(params.fasta_file_path).exists()) {
    exit 1, """
    Error: Missing or incorrect path for --fasta_file_path parameter.
    Please ensure the genome FASTA file exists.
    """
}
if (!params.out_dir || !new File(params.out_dir).exists()) {
    exit 1, """
    Error: Missing or incorrect path for --out_dir parameter.
    Please ensure the output directory exists.
    """
}
// Check if the BAM directory contains any .bam files
if (!new File(params.bam_directory_path).list().any { it.endsWith('.bam') }) {
    exit 1, """
    Error: No .bam files found in the input BAM directory.
    Please verify the path and content for --bam_directory_path.
    """
}
// Check if AT and GC anchors are within expected numeric ranges
if (!(params.at_anchor > 0 && params.at_anchor < 50)) {
    exit 1, """
    Error: at_anchor should be > 0 and < 50.
    Please check your input.
    """
}
if (!(params.gc_anchor > 50 && params.gc_anchor < 100)) {
    exit 1, """
    Error: gc_anchor should be > 50 and < 100.
    Please check your input.
    """
}
// Check failure and warning thresholds
if (!(params.failure_fold_change > 0)) {
    exit 1, """
    Error: failure_fold_change should be > 0.
    Please check your input.
    """
}
if (!(params.warning_fold_change > 0)) {
    exit 1, """
    Error: warning_fold_change should be > 0.
    Please check your input.
    """
}
if (params.warning_fold_change > params.failure_fold_change) {
    exit 1, """
    Error: failure_fold_change should be higher than warning_fold_change.
    Please check your input.
    """
}
if (!(params.failure_at > 0 && params.failure_gc > 0)) {
    exit 1, """
    Error: Both failure_at and failure_gc should be > 0.
    Please check your input.
    """
}
// Check boolean parameters
if (!(params.draw_trend.toString().toLowerCase() in ['true', '1', 't', 'false', '0', 'f'])) {
    exit 1, """
    Error: --draw_trend should be boolean: true or false.
    Please check your input.
    """
}
if (!(params.show_sample_names.toString().toLowerCase() in ['true', '1', 't', 'false', '0', 'f'])) {
    exit 1, """
    Error: --show_sample_names should be boolean: true or false.
    Please check your input.
    """
}

reference_fasta = file(params.fasta_file_path)
probe_bed = file(params.bed_file_path)

// Extract library names and pair with files
bam_channel = Channel
                .fromPath("${params.bam_directory_path}/*.bam")
                .map { file -> 
                    def libraryName = file.getName().replace(".bam", "")
                    return tuple(libraryName, file) 
                }

process bedtools_intersect {
    /*
     * Run bedtools intersect
     */
    input:
    tuple val(library_id), path(input_bam)

    output:
    tuple val(library_id), path ("${library_id}_intersected.bam")

    script:
    """
    bedtools intersect -a $input_bam -b $probe_bed > ${library_id}_intersected.bam
    """
}

process bedtools_coverage {
    maxForks 4
    /*
     * Run bedtools coverage
     */

    input:
    tuple val(library_id), path(intersected_bam)

    output:
    path("${library_id}_intersected_coverage.bed")

    script:
    """
    cut -f1-3 $probe_bed >> only_required_columns.bed
    bedtools coverage -d -a only_required_columns.bed -b $intersected_bam > ${library_id}_intersected_coverage.bed
    """
}

process bedtools_getfasta {
    /*
     * Run bedtools getfasta
     */
    input:
    path probe_bed
    path reference_fasta

    output:
    path "extracted_sequences.fasta"

    script:
    """
    bedtools getfasta -fi $reference_fasta -bed $probe_bed -fo extracted_sequences.fasta
    """
}

process calculate_gc_content {
    /*
     * Run helper bash script to create GC content.
     */
    input:
    path input_extracted_sequences

    output:
    path "extracted_sequences_GC_content.txt"

    script:
    """
    calculate_gc_content.sh $input_extracted_sequences extracted_sequences_GC_content.txt
    """
}

process create_soft_links {
    // Inputs: a list of file paths
    input:
    path bed_files

    // Output: path of a temporary directory containing soft links
    output:
    path "temporary_folder"

    script:
    """
    # Create a temporary directory
    mkdir -p "temporary_folder"
    
    # Create soft links for each file in the directory
    for file in ${bed_files}
    do
        echo \$file
        ln -s \$(readlink -f \$file) "temporary_folder/"
    done
    """
}

process generate_gc_bias {
    /*
     * Run panelGC_main.R executable
     */
    input:
    path probe_bed
    path intersected_coverage_dir
    path gc_content_summary
    path out_dir
    val at_anchor
    val gc_anchor
    val failure_fold_change
    val warning_fold_change
    val failure_at
    val failure_gc
    val draw_trend
    val show_sample_names

    script:
    """
    panelGC_main.R --probe_bed_file $probe_bed --bam_coverage_directory $intersected_coverage_dir \
    --reference_gc_content_file $gc_content_summary --outdir $out_dir \
    --at_anchor $at_anchor --gc_anchor $gc_anchor \
    --failure_fold_change $failure_fold_change --warning_fold_change $warning_fold_change \
    --failure_at $failure_at --failure_gc $failure_gc \
    --draw_trend $draw_trend --show_sample_names $show_sample_names
    """
}

workflow {
    coverage_files = bam_channel
    			| bedtools_intersect
    			| bedtools_coverage
    
   
    gc_content_summary = bedtools_getfasta(probe_bed, reference_fasta) 
			| calculate_gc_content

    
    generate_gc_bias(
		probe_bed,
		create_soft_links(coverage_files.collect()),
		gc_content_summary,
		file(params.out_dir),
        params.at_anchor,
        params.gc_anchor,
        params.failure_fold_change,
        params.warning_fold_change,
        params.failure_at,
        params.failure_gc,
        params.draw_trend,
        params.show_sample_names
   )
}

workflow.onError {
    error_msg = """
    Error: Pipeline execution stopped.
    Please ensure the right genome (--fasta_file_path) and probes/genomic bins (--bed_file_path) are used.
    Please ensure the right alignment bams (--bam_directory_path) are supplied.
    If the error persists, increase --at_anchor and/or decrease --gc_anchor, as your dataset may have a narrower GC content distribution.
    """
    println error_msg
}