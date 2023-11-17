#!/usr/bin/env nextflow

nextflow.enable.dsl=2
calculate_gc_content_script = "$PWD/calculate_gc_content.sh"

params.bam_directory_path = ''
params.bed_file_path = ''
params.fasta_file_path = ''

// Check if the paths are not empty and the files exist
if (!params.bam_directory_path || !new File(params.bam_directory_path).exists()) {
    exit 1, "Missing or incorrect path for --bam_directory_path parameter!"
}
if (!params.bed_file_path || !new File(params.bed_file_path).exists()) {
    exit 1, "Missing or incorrect path for --bed_file_path parameter!"
}
if (!params.fasta_file_path || !new File(params.fasta_file_path).exists()) {
    exit 1, "Missing or incorrect path for --fasta_file_path parameter!"
}

// Check if the BAM directory contains any .bam files
if (!new File(params.bam_directory_path).list().any { it.endsWith('.bam') }) {
    exit 1, "The specified BAM directory does not contain any .bam files!"
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
    /*
     * Run bedtools coverage
     */
    input:
    tuple val(library_id), path(intersected_bam)

    output:
    tuple val(library_id), path("${library_id}_intersected_coverage.bed")

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
    bedtools getfasta -fi $reference_fasta -bed $probe_bed > extracted_sequences.fasta
    """
}

process generate_gc_content_summary {
    /*
     * Run helper bash script to create GC content.
     */
    input:
    path input_extracted_sequences

    output:
    path "extracted_sequences_GC_content.txt"

    script:
    """
    $calculate_gc_content_script $input_extracted_sequences extracted_sequences_GC_content.txt
    """
}

workflow {
    bam_channel
    | bedtools_intersect
    | bedtools_coverage

    bedtools_getfasta(probe_bed, reference_fasta) | generate_gc_content_summary
}
