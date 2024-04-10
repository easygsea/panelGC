#! /usr/bin/env Rscript
message("Loading packages and functions...")

suppressPackageStartupMessages({
  library(argparser)
  library(data.table)
  library(GenomicRanges)
  library(readxl)
  library(rtracklayer)
  library(tidyverse)
})

# Suppress summarise info.
options(dplyr.summarise.inform = FALSE)

# Global variables ----
GC_UPPER_ANCHOR <- 75
GC_LOWER_ANCHOR <- 25
RELATIVE_FOLD_CHANGE_FAILURE_CUTOFF <- 2
RELATIVE_FOLD_CHANGE_WARNING_CUTOFF <- 1.5

BIAS_COLORS <- c(
  "GC biased" = "blue",
  "AT biased" = "red",
  "GC bias warning" = "cornflowerblue",
  "AT bias warning" = "coral",
  "no bias" = "darkgrey"
)
read_coverage_files_and_create_tibbles <- function(
    directory_path,
    file_pattern) {
  # List all files in the directory that match the given pattern
  file_paths <- list.files(
    path = directory_path,
    pattern = file_pattern,
    full.names = TRUE
  )

  # Initialize an empty list to store tibbles
  tibbles_list <- list()

  # Read each file, create a tibble, and assign it a name based on the file name
  for (file_path in file_paths) {
    # Extract the base name of the file (without path and extension)
    file_name <- tools::file_path_sans_ext(basename(file_path))

    # Read the file into a tibble
    tibble_data <- read_tsv(file_path,
      col_names = c(
        "chromosome",
        "probe_start_pos",
        "probe_end_pos",
        "offset",
        "depth"
      ),
      col_types = cols(
        chromosome = col_character(),
        probe_start_pos = col_double(),
        probe_end_pos = col_double(),
        offset = col_double(),
        depth = col_double()
      )
    )

    # Assign the tibble to the list with the name as the key
    tibbles_list[[file_name]] <- tibble_data
  }

  return(tibbles_list)
}

read_gc_content_file <- function(file_path) {
  # Read the TSV file
  data <- read_tsv(
    file_path,
    col_types = cols(
      region = col_character(),
      GC = col_double()
    )
  )

  # Process the data
  reformatted_data <- data %>%
    separate(region, into = c("chromosome", "positions"), sep = "[:]") %>%
    separate(positions, into = c("start_pos", "end_pos"), sep = "-") %>%
    mutate(
      start_pos = as.integer(start_pos),
      end_pos = as.integer(end_pos)
    )

  return(reformatted_data)
}

read_bed_file <- function(file_path) {
  message("Begin reading bed file content.")
  # Supply your BED file path with --bed_file_path
  bed_data <- suppressMessages(
    read_tsv(file_path, col_names = FALSE)
  )

  # Check the number of columns and process accordingly
  if (ncol(bed_data) >= 4) {
    colnames(bed_data)[4] <- "probe_name"
  } else {
    bed_data$probe_name <- seq(1000, 1000 + nrow(bed_data) - 1)
  }

  # Ensuring the first three columns are named correctly
  colnames(bed_data)[1:3] <- c("chromosome", "probe_start_pos", "probe_end_pos")
  message("End reading bed content.")
  
  return(bed_data)
}
add_offset_to_start_position_and_drop_columns <- function(tibble) {
  tibble %>%
    # Add offset to start position. Subtract one as offset starts with 1.
    mutate(position = probe_start_pos + offset - 1) %>%
    # Select the required columns
    select(sample, probe_name, chromosome, position, depth)
}
calculate_gc_bias_regression <- function(bin_gc_summary) {
  gc_bias_regression <- bin_gc_summary %>%
    group_by(sample) %>%
    arrange(sample, gc_percentile) %>%
    nest() %>%
    mutate(loess_depth = purrr::map(data, function(x) {
      stats::loess(normalized_depth ~ gc_percentile, span = 0.75, data = x) %>%
        stats::predict(gc_percentile = unique(gc_percentile))
    })) %>%
    unnest(cols = c(data, loess_depth))
  return(gc_bias_regression)
}

get_gc_bias_regression_table <- function(
    raw_bam_readcount_intersected_probes,
    gc_content) {
  mean_sample_depths <- raw_bam_readcount_intersected_probes %>%
    select(sample, chromosome, position, depth) %>%
    distinct() %>%
    group_by(sample) %>%
    summarise(mean_depth = mean(depth))

  # Compute the median depth of each probe (or genomic bin) per sample.
  gc_summary <- raw_bam_readcount_intersected_probes %>%
    inner_join(gc_content,
      multiple = "all"
    ) %>%
    group_by(probe_name, GC, sample) %>%
    summarise(depth_median = median(depth))

  # Compute the median depth across probes (or genomic bins) with the same GC content.
  bin_gc_summary <- gc_summary %>%
    group_by(sample) %>%
    mutate(gc_percentile = floor(GC * 100) / 100) %>%
    group_by(sample, gc_percentile) %>%
    summarise(depth_median_median = median(depth_median))

  # Normalization by mean sample depth.
  bin_gc_summary <- left_join(
    bin_gc_summary, mean_sample_depths,
    by = "sample"
  ) %>%
    mutate(normalized_depth = log2(
      depth_median_median / mean_depth + 1
    ))

  gc_bias_regression <- calculate_gc_bias_regression(bin_gc_summary)
}
compute_bias_cutoff <- function(fold_change_cutoff, is_upper) {
  bias_factor <- if (is_upper) fold_change_cutoff else 1 / fold_change_cutoff
  return(log2(bias_factor + 1))
}

compute_cutoffs <- function(cutoff) {
  upper <- compute_bias_cutoff(cutoff, TRUE)
  lower <- compute_bias_cutoff(cutoff, FALSE)
  return(list(upper = upper, lower = lower))
}

get_gc_bias_classification_table <- function(
    gc_bias_regression,
    relative_bias_cutoff_failure_upper,
    relative_bias_cutoff_failure_lower,
    relative_bias_cutoff_warning_upper,
    relative_bias_cutoff_warning_lower,
    relative_bias_name, at_bias_name,
    gc_bias_name) {
  gc_bias_loess <- calculate_gc_bias_loess(
    gc_bias_regression,
    relative_bias_name,
    at_bias_name,
    gc_bias_name
  )
  gc_bias_classification <- classify_gc_bias(
    gc_bias_loess,
    relative_bias_cutoff_failure_upper,
    relative_bias_cutoff_failure_lower,
    relative_bias_cutoff_warning_upper,
    relative_bias_cutoff_warning_lower,
    relative_bias_name,
    at_bias_name,
    gc_bias_name
  )
  return(gc_bias_classification)
}

classify_gc_bias <- function(
    gc_bias_loess,
    relative_bias_cutoff_failure_upper,
    relative_bias_cutoff_failure_lower,
    relative_bias_cutoff_warning_upper,
    relative_bias_cutoff_warning_lower,
    relative_bias_name,
    at_bias_name,
    gc_bias_name) {
  gc_bias_classification <- gc_bias_loess %>%
    pivot_longer(
      cols = matches("^bias_.*[0-9]$"),
      names_to = "bias_category",
      values_to = "bias_score"
    ) %>%
    filter(bias_category == relative_bias_name) %>%
    mutate(
      bias_type = case_when(
        bias_score >= relative_bias_cutoff_failure_upper ~ "GC biased",
        bias_score <= relative_bias_cutoff_failure_lower ~ "AT biased",
        bias_score >= relative_bias_cutoff_warning_upper ~ "GC bias warning",
        bias_score <= relative_bias_cutoff_warning_lower ~ "AT bias warning",
        TRUE ~ "no bias"
      )
    ) %>%
    # Convert Bias_Score into fold change.
    mutate(
      bias_score_fold_change = 2 ** bias_score - 1, .after = "bias_score"
    )
  return(gc_bias_classification)
}

# TODO: think about a name change to justify including
# this function in this module.
calculate_gc_bias_loess <- function(
    gc_bias_regression,
    relative_bias_name,
    at_bias_name,
    gc_bias_name) {
  gc_bias_loess <- gc_bias_regression %>%
    filter(
      gc_percentile == GC_UPPER_ANCHOR / 100 | gc_percentile == GC_LOWER_ANCHOR / 100
    ) %>%
    select(sample, gc_percentile, loess_depth) %>%
    mutate(gc_percentile = gc_percentile * 100) %>%
    pivot_wider(
      names_from = gc_percentile, values_from = loess_depth,
      names_prefix = "bias_"
    ) %>%
    # Compute GC-to-AT relative bias score.
    mutate(!!relative_bias_name :=
      log2((2 ** get(gc_bias_name) - 1) /
        (2 ** get(at_bias_name) - 1) + 1))
  return(gc_bias_loess)
}
plot_gc_profiles <- function(gc_bias_regression, gc_bias_classification) {
  # Maximum of y-axis.
  y_max <- ceiling(max(pull(gc_bias_regression, normalized_depth)))
  # Generate sample GC profiles plot.
  p <- gc_bias_regression %>%
    left_join(
      select(
        gc_bias_classification,
        "sample", "bias_type"
      ),
      by = "sample", multiple = "all"
    ) %>%
    ggplot(aes(x = gc_percentile, y = normalized_depth, color = bias_type)) +
    geom_smooth(aes(group = sample),
      method = "loess", formula = y ~ x, se = FALSE, fullrange = TRUE
    ) +
    scale_y_continuous("LOESS Depth Per GC Percentile") +
    ggtitle("GC Content vs. Coverage by Sample") +
    coord_cartesian(xlim = c(0.0, 1.0), ylim = c(0, y_max)) +
    scale_x_continuous("GC", breaks = seq(0, 1.0, 0.2)) +
    scale_color_manual(name = "Bias Type", values = BIAS_COLORS) +
    theme_bw(base_size = 20) +
    theme(legend.position = "right")
  return(p)
}

main <- function(
    probe_bed_file,
    bam_coverage_directory,
    reference_gc_content_file,
    outdir) {
  message("In main function")
  # Set variable names.
  relative_bias_name <- str_glue("bias_{GC_LOWER_ANCHOR}vs{GC_UPPER_ANCHOR}")
  at_bias_name <- str_glue("bias_{GC_LOWER_ANCHOR}")
  gc_bias_name <- str_glue("bias_{GC_UPPER_ANCHOR}")

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
  ## Read coverage data for each sample.
  raw_coverage_tibbles <- read_coverage_files_and_create_tibbles(
    bam_coverage_directory,
    "_intersected_coverage\\.bed$"
  )
  all_libraries_raw_coverage <- raw_coverage_tibbles %>%
    imap(~ mutate(.x, sample = .y)) %>%
    bind_rows() %>%
    mutate(sample = sub("_intersected_coverage$", "", sample))
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
  # else compute the relative HQ depth for each GC bin per sample, perform
  # LOESS regression, and obtain predicted relative HQ depth for each GC bin.
  write.table(
    gc_bias_regression_table,
    file = outfile_gc_bias_regression,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cutoffs_failure <- compute_cutoffs(RELATIVE_FOLD_CHANGE_FAILURE_CUTOFF)
  cutoff_failure_upper <- cutoffs_failure$upper
  cutoff_failure_lower <- cutoffs_failure$lower

  cutoffs_warning <- compute_cutoffs(RELATIVE_FOLD_CHANGE_WARNING_CUTOFF)
  cutoff_warning_upper <- cutoffs_warning$upper
  cutoff_warning_lower <- cutoffs_warning$lower

  gc_bias_classification_table <- get_gc_bias_classification_table(
    gc_bias_regression_table,
    cutoff_failure_upper,
    cutoff_failure_lower,
    cutoff_warning_upper,
    cutoff_warning_lower,
    relative_bias_name,
    at_bias_name,
    gc_bias_name
  )

  # The file recording LOESS regression results.
  outfile_gc_bias_classification <- file.path(
    outdir,
    "gc_bias_loess_classification.tsv"
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
      "R script for calculating and visualizing GC biases. Required inputs are",
      "bam files in a directory, probes BED file and reference sequence FASTA.",
      "Usage = generate_gc_biases.R --probe_bed_file <probe_bed_file>",
      "--bam_coverage_directory <bam_coverage_directory>",
      "--reference_gc_content_file <reference_gc_content_file>",
      "--outdir <outdir>"
    )
  )
  parser <- add_argument(parser, "--probe_bed_file", help = "Probes bed file")
  parser <- add_argument(
    parser,
    "--bam_coverage_directory",
    help = "Directory of coverage files for each sample. Expects
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