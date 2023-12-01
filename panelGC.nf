#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.bam_directory_path = ''
params.bed_file_path = ''
params.fasta_file_path = ''
params.out_dir = ''

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
if (!params.out_dir || !new File(params.out_dir).exists()) {
    exit 1, "Missing or incorrect path for --out_dir parameter!"
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
    calculate_gc_content.sh $input_extracted_sequences extracted_sequences_GC_content.txt
    """
}

process create_soft_links {
    // Inputs: a list of file paths
    input:
    path bed_files

    // Output: path of the new directory containing soft links
    output:
    path "temp_folder"

    script:
    """
    # Create a unique directory
    # Not sure if -p is essential
    mkdir -p "temp_folder"
    
    # Create soft links for each file in the directory
    for file in ${bed_files}
    do
        echo \$file
        ln -s \$(readlink -f \$file) "temp_folder/"
    done
    """
}

process generate_gc_bias {
    /*
     * Run generate_gc_bias executable
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
    source /projects/clingenetics/mgoktas_dev/ccgi_standalone_utilities/R_utilities/R_profile.sh 
    panelGC_main.R --probe_bed_file $probe_bed --bam_coverage_directory $intersected_coverage_dir --reference_gc_content_file $gc_content_summary --outdir $out_dir
    """
}

workflow {
    coverage_files = bam_channel
    			| bedtools_intersect
    			| bedtools_coverage
    
   
    gc_content_summary = bedtools_getfasta(probe_bed, reference_fasta) 
			| generate_gc_content_summary

    
    generate_gc_bias(
		probe_bed,
		create_soft_links(coverage_files.collect()),
		gc_content_summary,
		params.out_dir
   )

}
