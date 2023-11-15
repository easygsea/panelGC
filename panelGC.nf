#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.bam_path = ''
params.bed = ''
params.fasta = ''

if (!params.bam_path || !params.bed || !params.fasta) {
    exit 1, "Missing input files: --bam, --bed, and --fasta parameters are required!"
}

input_fasta = file(params.fasta)
input_bed = file(params.bed)

calculate_gc_content_script = "$PWD/calculate_gc_content.sh"

// Extract library names and pair with files
bam_channel = Channel
                .fromPath(params.bam_path + '/*.bam')
                .map { file -> 
                    def libraryName = file.getName().replace(".bam", "")
                    return tuple(libraryName, file) 
                }

bam_channel.view{ tuple -> 
    println "Library: ${tuple[0]}, File: ${tuple[1]}"
}

process IntersectBedtools {
    /*
     * Run bedtools intersect
     */
    input:
    tuple val(library_id), path(input_bam)

    output:
    tuple val(library_id), path ("${library_id}_intersected.bam")

    script:
    """
    bedtools intersect -a $input_bam -b $input_bed > ${library_id}_intersected.bam
    """
}

process CoverageBedtools {
    /*
     * Run bedtools coverage
     */
    input:
    tuple val(library_id), path(intersected_bam)

    output:
    tuple val(library_id), path("${library_id}_intersected_coverage.bed")

    script:
    """
    bedtools coverage -a $input_bed -b $intersected_bam > ${library_id}_intersected_coverage.bed
    """
}

process GetFastaBedtools {
    /*
     * Run bedtools getfasta
     */
    input:
    path input_bed
    path input_fasta

    output:
    path "extracted_sequences.fasta"

    script:
    """
    bedtools getfasta -fi $input_fasta -bed $input_bed > extracted_sequences.fasta
    """
}

process GENERATE_GC_CONTENT{
    /*
     * Run helper bash script to create GC content.
     */
    input:
    path input_extracted_sequences

    output:
    path "extracted_sequences_GC_content.txt"

    script:
    """
    bash $calculate_gc_content_script $input_extracted_sequences extracted_sequences_GC_content.txt
    """ 

}

workflow {
    intersect_receiver_ch = IntersectBedtools(bam_channel)
    coverage_bedtools_receiver_ch = CoverageBedtools(intersect_receiver_ch)
    get_fasta_bedtools_receiver_channel = GetFastaBedtools(input_bed, input_fasta)
    generate_gc_content_receiver_channel = GENERATE_GC_CONTENT(get_fasta_bedtools_receiver_channel)
}
