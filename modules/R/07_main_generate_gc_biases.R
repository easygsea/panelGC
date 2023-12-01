#! /usr/bin/env Rscript

source("01_load_libraries_and_constants.R")
source("02_load_data.R")
source("03_data_manipulators.R")
source("04_gc_bias_regression.R")
source("05_gc_bias_classification.R")
source("06_gc_bias_plot.R")

main <- function(
  probe_bed_file,
  bam_coverage_directory,
  reference_gc_content_file,
  outdir) {
  message("In main function")
  # Set variable names.
  relative_bias_name <- str_glue("bias_{GC.LOWER.ANCHOR}vs{GC.UPPER.ANCHOR}")
  at_bias_name <- str_glue("bias_{GC.LOWER.ANCHOR}")
  gc_bias_name <- str_glue("bias_{GC.UPPER.ANCHOR}")

  ## ------ Read files.
  ## Read probes.
  probes <- read_bed_file(probe_bed_file)

  ## Read GC content
  reference_gc_content <- read_gc_content_file(reference_gc_content_file)
  ## Add probe name/identifier to reference_gc_content.
  reference_gc_content <- reference_gc_content %>% left_join(
    probes,
    by = c(
      "chromosome" = "chromosome",
      "start_pos" = "probe_start_pos",
      "end_pos" = "probe_end_pos"
    )
  )
  ## Read coverage data for each library.
  raw_coverage_tibbles <- read_coverage_files_and_create_tibbles(
    bam_coverage_directory,
    "_intersected_coverage\\.bed$"
  )
  all_libraries_raw_coverage <- raw_coverage_tibbles %>%
    imap(~ mutate(.x, library = .y)) %>%
    bind_rows()
  ## Add probe name/identifier to reference_gc_content.
  all_libraries_raw_coverage <- all_libraries_raw_coverage %>% left_join(
    probes,
    by = c(
      "chromosome" = "chromosome",
      "probe_start_pos" = "probe_start_pos",
      "probe_end_pos" = "probe_end_pos"
    )
  )
  all_libraries_raw_coverage <- add_offset_to_start_position_and_drop_columns(
    all_libraries_raw_coverage
  )

  gc_bias_regression_table <- get_gc_bias_regression_table(
    all_libraries_raw_coverage,
    reference_gc_content
  )
  # The file recording LOESS regression results.
  outfile_gc_bias_regression <- file.path(
    outdir, "gc_bias_loess_regression.tsv"
  )
  # Check if LOESS regression table exists and overwrite disabled. If yes, load;
  # else compute the relative HQ depth for each GC bin per library, perform
  # LOESS regression, and obtain predicted relative HQ depth for each GC bin.
  write.table(
    gc_bias_regression_table,
    file = outfile_gc_bias_regression,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cutoffs_failure <- compute_cutoffs(RELATIVE.FOLD.CHANGE.FAILURE.CUTOFF)
  cutoff_failure_upper <- cutoffs_failure$upper
  cutoff_failure_lower <- cutoffs_failure$lower

  cutoffs_warning <- compute_cutoffs(RELATIVE.FOLD.CHANGE.WARNING.CUTOFF)
  cutoff_warning_upper <- cutoffs_warning$upper
  cutoff_warning_lower <- cutoffs_warning$lower

  gc_bias_classification_table <- get_gc_bias_classification_table(
    gc_bias_regression_table,
    cutoff_failure_upper, cutoff_failure_lower,
    cutoff_warning_upper, cutoff_warning_lower,
    relative_bias_name, at_bias_name, gc_bias_name
  )

  # The file recording LOESS regression results.
  outfile_gc_bias_classification <- file.path(
    outdir, "gc_bias_loess_classification.tsv"
  )

  write.table(
    gc_bias_classification_table,
    file = outfile_gc_bias_classification,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  print(gc_bias_regression_table)
  print(gc_bias_classification_table)
  p <- plot_gc_profiles(gc_bias_regression_table, gc_bias_classification_table)
  ggsave(
    file.path(outdir, "gc_bias_curves.png"),
    plot = p,
    height = 7.5,
    width = 10,
    units = "in"
  )
}

# * Functions: args ----
parse_args_function <- function() {
  parser <- arg_parser(
    paste(
      "R script for calculating and visualizing library GC biases for each",
      "library in the input directory. This program requires bam fileof each",
      "library in a directory, probes BED file and reference sequence FASTA",
      "file. Usage = generate_gc_biases.R --probe_bed_file <probe_bed_file>",
      "--bam_coverage_directory <bam_coverage_directory>",
      "--reference_gc_content_file <reference_gc_content_file>",
      "--outdir <outdir>"
    )
  )
  parser <- add_argument(parser, "--probe_bed_file", help = "Probes bed file")
  parser <- add_argument(
    parser,
    "--bam_coverage_directory",
    help = "Directory of coverage files for each library. Expects
	'_intersected_coverage.bed'	suffix on each file."
  )
  parser <- add_argument(
    parser,
    "--reference_gc_content_file",
    help = "Probe sequence GC content tab seperated file."
  )
  parser <- add_argument(
    parser, "--outdir",
    help = "Pool analysis directory"
  )
  argv <- parse_args(parser)
  return(argv)
}

if (!interactive()) {
  args <- parse_args_function()
  probe_bed_file <- pluck(args, "probe_bed_file")
  bam_coverage_directory <- pluck(args, "bam_coverage_directory")
  reference_gc_content_file <- pluck(args, "reference_gc_content_file")
  outdir <- pluck(args, "outdir")
  main(
    probe_bed_file,
    bam_coverage_directory,
    reference_gc_content_file,
    outdir
  )
}
