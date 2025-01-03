#! /usr/bin/env Rscript
suppressPackageStartupMessages({
  library(argparser)
  library(GenomicRanges)
  library(rtracklayer)
  library(tidyverse)
})

# Suppress summarise info.
options(dplyr.summarise.inform = FALSE)

# Global variables ----
BIAS_COLORS <- c(
  "GC biased" = "blue",
  "AT biased" = "red",
  "GC bias warning" = "cornflowerblue",
  "AT bias warning" = "coral",
  "no bias" = "darkgrey"
)

read_sample_labels_csv <- function(file_path) {
  sample_labels <- read_csv(
    file_path,
    col_names = TRUE,
    col_types = cols(
      .default = col_character()
    )
  )
  if(colnames(sample_labels)[1] != "sample"){
    stop("First column must be named 'sample' in sample_labels_csv.")
  }
  return(sample_labels)
}

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
    cutoff_failure_at,
    cutoff_failure_gc,
    relative_bias_name,
    at_bias_name,
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
    cutoff_failure_at,
    cutoff_failure_gc,
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
    cutoff_failure_at,
    cutoff_failure_gc,
    relative_bias_name,
    at_bias_name,
    gc_bias_name) {
  relative_bias_symbol <- as.symbol(relative_bias_name)
  at_bias_symbol <- as.symbol(at_bias_name)
  gc_bias_symbol <- as.symbol(gc_bias_name)
  relative_failture_threshold_at_name <- str_glue("{relative_bias_name}_failure_threshold_at")
  relative_failture_threshold_gc_name <- str_glue("{relative_bias_name}_failure_threshold_gc")
  relative_warning_threshold_at_name <- str_glue("{relative_bias_name}_warning_threshold_at")
  relative_warning_threshold_gc_name <- str_glue("{relative_bias_name}_warning_threshold_gc")
  at_threshold_name <- str_glue("{at_bias_name}_failure_threshold")
  gc_threshold_name <- str_glue("{gc_bias_name}_failure_threshold")
  gc_bias_classification <- gc_bias_loess %>%
    select(
      c(
        "sample",
        all_of(relative_bias_name),
        all_of(at_bias_name),
        all_of(gc_bias_name)
      )
    ) %>%
    # Convert bias score into fold change.
    mutate(
      across(
        starts_with("bias_"), ~ 2 ** . - 1,
        .names = "{.col}_fold_change"
      )
    ) %>%
    mutate(
      !!relative_failture_threshold_at_name := relative_bias_cutoff_failure_lower,
      !!relative_failture_threshold_gc_name := relative_bias_cutoff_failure_upper,
      !!relative_warning_threshold_at_name := relative_bias_cutoff_warning_lower,
      !!relative_warning_threshold_gc_name := relative_bias_cutoff_warning_upper,
      !!at_threshold_name := cutoff_failure_at,
      !!gc_threshold_name := cutoff_failure_gc,
      bias_type = case_when(
        !!relative_bias_symbol >= relative_bias_cutoff_failure_upper ~ "GC biased",
        !!relative_bias_symbol <= relative_bias_cutoff_failure_lower ~ "AT biased",
        !!gc_bias_symbol >= cutoff_failure_gc ~ "GC biased",
        !!at_bias_symbol >= cutoff_failure_at ~ "AT biased",
        !!relative_bias_symbol >= relative_bias_cutoff_warning_upper ~ "GC bias warning",
        !!relative_bias_symbol <= relative_bias_cutoff_warning_lower ~ "AT bias warning",
        TRUE ~ "no bias"
      )
    ) %>%
    select(
      c(
        sample,
        starts_with(relative_bias_name),
        starts_with(at_bias_name),
        starts_with(gc_bias_name),
        bias_type
      )
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
      gc_percentile == GC_ANCHOR / 100 | gc_percentile == AT_ANCHOR / 100
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

plot_gc_profiles <- function(gc_bias_regression, gc_bias_classification, sample_labels) {
  label <- NULL
  gc_bias_regression_w_labels <- gc_bias_regression
  # Check if sample_labels file is provided by the user
  if (!is.null(sample_labels)) {
    # Second column name is the label
    label <- colnames(sample_labels)[2]
    gc_bias_regression_w_labels <- gc_bias_regression %>%
      left_join(
        sample_labels,
        by = "sample"
      )
  }

  # Maximum of y-axis.
  y_max <- ceiling(max(pull(gc_bias_regression, normalized_depth)))
  # Generate sample GC profiles plot.
  base_plot <- gc_bias_regression_w_labels %>%
    left_join(
      select(
        gc_bias_classification,
        "sample", "bias_type"
      ),
      by = "sample", multiple = "all"
    ) %>%
    ggplot(aes(x = gc_percentile, y = normalized_depth, color = bias_type))

  # Conditionally add the geom_smooth layer
  if (!is.null(label)) {
    geom_smooth_layer <- geom_smooth(aes(group = sample, linetype = .data[[label]]),
                                     method = "loess", formula = y ~ x, se = FALSE, fullrange = TRUE)
  } else {
    geom_smooth_layer <- geom_smooth(aes(group = sample),
                                     method = "loess", formula = y ~ x, se = FALSE, fullrange = TRUE)
  }

  # Add the rest of the layers
  p <- base_plot +
    geom_smooth_layer +
    scale_y_continuous("LOESS Depth Per GC Percentile") +
    ggtitle("GC Content vs. Coverage by Sample") +
    coord_cartesian(xlim = c(0.0, 1.0), ylim = c(0, y_max)) +
    scale_x_continuous("GC Content", breaks = seq(0, 1.0, 0.2)) +
    scale_color_manual(name = "Bias Type", values = BIAS_COLORS) +
    theme_bw(base_size = 20) +
    theme(legend.position = "right")

  return(p)
}

main <- function(
    probe_bed_file,
    bam_coverage_directory,
    reference_gc_content_file,
    sample_labels_csv,
    outdir) {
  # Set variable names.
  relative_bias_name <- str_glue("bias_{AT_ANCHOR}vs{GC_ANCHOR}")
  at_bias_name <- str_glue("bias_{AT_ANCHOR}")
  gc_bias_name <- str_glue("bias_{GC_ANCHOR}")

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

  ## Read sample types if provided
  sample_labels <- if (sample_labels_csv != "NA") {
    read_sample_labels_csv(sample_labels_csv)
  } else {
    NULL
  }

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
  write.table(
    gc_bias_regression_table,
    file = outfile_gc_bias_regression,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cutoffs_failure <- compute_cutoffs(FAILURE_FOLD_CHANGE)
  cutoff_failure_upper <- cutoffs_failure$upper
  cutoff_failure_lower <- cutoffs_failure$lower

  cutoffs_warning <- compute_cutoffs(WARNING_FOLD_CHANGE)
  cutoff_warning_upper <- cutoffs_warning$upper
  cutoff_warning_lower <- cutoffs_warning$lower
  
  cutoff_failure_at <- compute_bias_cutoff(FAILURE_AT, TRUE)
  cutoff_failure_gc <- compute_bias_cutoff(FAILURE_GC, TRUE)

  gc_bias_classification_table <- get_gc_bias_classification_table(
    gc_bias_regression_table,
    cutoff_failure_upper,
    cutoff_failure_lower,
    cutoff_warning_upper,
    cutoff_warning_lower,
    cutoff_failure_at,
    cutoff_failure_gc,
    relative_bias_name,
    at_bias_name,
    gc_bias_name
  )

  # The file recording bias classification results.
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

  p <- plot_gc_profiles(gc_bias_regression_table, gc_bias_classification_table, sample_labels)
  ggsave(
    file.path(outdir, "gc_bias_profile.png"),
    plot = p,
    height = 7.5,
    width = 10,
    units = "in"
  )
  
  # Draw trend visualization.
  if(DRAW_TREND) {
    source(file.path(BIN_FOLDER, "modules", "plot_trend_function.R"))
    p <- plot_gc_trend(gc_bias_classification_table, SHOW_SAMPLES)
    p_width <- ifelse(SHOW_SAMPLES, 9, 7)
    ggsave(
      file.path(outdir, "gc_bias_trend.png"),
      plot = p,
      height = 2 + nrow(gc_bias_classification_table) * 0.25,
      width = p_width,
      units = "in"
    )
  }
}

# * Functions: args ----
find_here <- function() {
  raw_args <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", raw_args[grep("--file=", raw_args)])
  dirname(normalizePath(script_path))
}

parse_args_function <- function() {
  parser <- arg_parser(
    paste(
      "R script for calculating and visualizing GC biases. Required inputs are",
      "bam files in a directory, probes BED file and reference sequence FASTA.",
      "Usage = generate_gc_biases.R --probe_bed_file <probe_bed_file>",
      "--bam_coverage_directory <bam_coverage_directory>",
      "--reference_gc_content_file <reference_gc_content_file>",
      "--sample_labels_csv <sample_labels_csv>",
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
    parser,
    "--sample_labels_csv",
    help = "Sample labels csv for grouping by line type in 'gc_bias_profile.png', optional",
    default = "NA"
  )
  parser <- add_argument(
    parser,
    "--outdir",
    help = "Output directory."
  )
  parser <- add_argument(
    parser,
    "--at_anchor",
    help = "GC percentile anchor for detecting AT bias."
  )
  parser <- add_argument(
    parser,
    "--gc_anchor",
    help = "GC percentile anchor for detecting GC bias."
  )
  parser <- add_argument(
    parser,
    "--failure_fold_change",
    help = "Relative coverage fold change failure threshold."
  )
  parser <- add_argument(
    parser,
    "--warning_fold_change",
    help = "Relative coverage fold change warning threshold."
  )
  parser <- add_argument(
    parser,
    "--failure_at",
    help = "Coverage fold change failure threshold at the AT anchor."
  )
  parser <- add_argument(
    parser,
    "--failure_gc",
    help = "Coverage fold change failure threshold at the GC anchor."
  )
  parser <- add_argument(
    parser,
    "--draw_trend",
    help = "Generate trend visualization."
  )
  parser <- add_argument(
    parser,
    "--show_sample_names",
    help = "Show sample names in trend visualization"
  )
  argv <- parse_args(parser)
  return(argv)
}

if (!interactive()) {
  args <- parse_args_function()
  probe_bed_file <- pluck(args, "probe_bed_file")
  bam_coverage_directory <- pluck(args, "bam_coverage_directory")
  reference_gc_content_file <- pluck(args, "reference_gc_content_file")
  sample_labels_csv <- pluck(args, "sample_labels_csv")
  outdir <- pluck(args, "outdir")
  AT_ANCHOR <<- floor(as.numeric(pluck(args, "at_anchor")))
  GC_ANCHOR <<- floor(as.numeric(pluck(args, "gc_anchor")))
  FAILURE_FOLD_CHANGE <<- as.numeric(pluck(args, "failure_fold_change"))
  WARNING_FOLD_CHANGE <<- as.numeric(pluck(args, "warning_fold_change"))
  FAILURE_AT <<- as.numeric(pluck(args, "failure_at"))
  FAILURE_GC <<- as.numeric(pluck(args, "failure_gc"))
  DRAW_TREND <<- as.logical(pluck(args, "draw_trend"))
  SHOW_SAMPLES <<- as.logical(pluck(args, "show_sample_names"))
  BIN_FOLDER <<- find_here()
  main(
    probe_bed_file,
    bam_coverage_directory,
    reference_gc_content_file,
    sample_labels_csv,
    outdir
  )
}
