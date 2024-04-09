#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.bam_directory_path = ''
params.bed_file_path = ''
params.fasta_file_path = ''
params.out_dir = ''

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
    bedtools coverage -d -abam $intersected_bam -b only_required_columns.bed > ${library_id}_intersected_coverage.bed
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

    output:
    path "${out_dir}/gc_bias_curves.png"

    script:
    """
    panelGC_main.R --probe_bed_file $probe_bed --bam_coverage_directory $intersected_coverage_dir \
    --reference_gc_content_file $gc_content_summary --outdir $out_dir
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
		file(params.out_dir)
   )

}