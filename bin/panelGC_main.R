#! /usr/bin/env Rscript
suppressPackageStartupMessages({
  library(argparser)
  library(data.table)
  library(tidyverse)
})

# Global variables ----
BIAS_COLORS <- c(
  "GC biased" = "blue",
  "AT biased" = "red",
  "GC bias warning" = "cornflowerblue",
  "AT bias warning" = "coral",
  "no bias" = "darkgrey"
)

read_sample_labels_csv <- function(file_path) {
  sample_labels <- fread(
    file_path,
    colClasses = list(
      character = names(fread(file_path, nrows = 0))
    )
  )
  if(colnames(sample_labels)[1] != "sample"){
    stop("First column must be named 'sample' in sample_labels_csv.")
  }
  return(sample_labels)
}

read_coverage_files_and_create_data_tables <- function(
    directory_path,
    file_pattern) {
  # List all files in the directory that match the given pattern
  file_paths <- list.files(
    path = directory_path,
    pattern = file_pattern,
    full.names = TRUE
  )

  # Initialize an empty list to store tibbles
  data_tables_list <- list()

  # Read each file, create a tibble, and assign it a name based on the file name
  for (file_path in file_paths) {
    # Extract the base name of the file (without path and extension)
    file_name <- tools::file_path_sans_ext(basename(file_path))

    # Read the file using data.table
    data_table <- fread(file_path,
      col.names = c(
        "chromosome",
        "region_start_pos", 
        "region_end_pos",
        "window_position",
        "depth"
      ),
      colClasses = list(
        character = "chromosome",
        numeric = c("region_start_pos", "region_end_pos", "window_position", "depth")
      )
    )[, region := sprintf("%s:%d-%d", chromosome, region_start_pos, region_end_pos)]

    # Assign the data table to the list with the name as the key
    data_tables_list[[file_name]] <- data_table
  }

  return(data_tables_list)
}

read_gc_content_file <- function(file_path) {
  # Read the TSV file using data.table
  data <- fread(
    file_path,
    colClasses = list(
      character = c("region", "region_name"),
      numeric = c("window_start", "window_end", "GC")
    )
  )
  return(data)
}

calculate_gc_bias_regression <- function(bin_gc_summary) {
  # Calculate LOESS regression for each sample
  gc_bias_regression <- bin_gc_summary[!is.na(normalized_depth), {
    setorder(.SD, gc_percentile)
    loess_model <- stats::loess(normalized_depth ~ gc_percentile,
                               span = 0.75, 
                               data = as.data.frame(.SD))

    # Create sequence of all possible GC percentiles from min to max.
    all_gc_percentiles <- seq(min(.SD$gc_percentile), max(.SD$gc_percentile), by = 0.01)

    # Predict for all percentiles.
    loess_depth <- stats::predict(loess_model,
                                 newdata = data.frame(gc_percentile = all_gc_percentiles))

    # Create data table with all percentiles and predictions.
    result <- data.table(
      gc_percentile = all_gc_percentiles,
      mean_depth = .SD$mean_depth,
      loess_depth = loess_depth
    )

    # Join with original data and clean up mean_depth columns
    result <- merge(result, .SD, by = "gc_percentile", all.x = TRUE)[
      , `:=`(mean_depth = mean_depth.x, mean_depth.x = NULL, mean_depth.y = NULL)
    ]
  }, by = sample]

  return(gc_bias_regression)
}

get_gc_bias_regression_table <- function(
    raw_bam_readcount_intersected_probes,
    gc_content) {
  # Calculate mean depths.
  mean_sample_depths <- raw_bam_readcount_intersected_probes[
    , position := region_start_pos + window_position - 1
  ][
    , .(sample, chromosome, position, depth)
  ][
    , .SD[1], by = .(sample, chromosome, position)
  ][
    , .(mean_depth = mean(depth)), by = sample
  ]

  # Expand GC content windows and join with read depth data.
  gc_summary <- gc_content[
    , .(window_position = seq.int(window_start, window_end)),
    by = .(region, region_name, window_start, window_end, GC)
  ][
    raw_bam_readcount_intersected_probes,
    on = c("region", "window_position"),
    mult = "all",
    nomatch = 0,  # Skip rows with no matches
    allow.cartesian = TRUE  # Allow cartesian join to handle multiple matches
  ][
    , .(depth_median = median(depth, na.rm = TRUE)),
    by = .(region, region_name, window_start, window_end, GC, sample)
  ]

  # Calculate GC percentiles and median depths.
  bin_gc_summary <- gc_summary[
    , .(depth_median_median = median(depth_median, na.rm = TRUE)),
    by = .(sample, gc_percentile = floor(GC * 100) / 100)
  ]

  # Join with mean depths and calculate normalized depth
  bin_gc_summary <- bin_gc_summary[
    mean_sample_depths,
    on = "sample"
  ][
    , normalized_depth := log2(depth_median_median / mean_depth + 1)
  ]

  gc_bias_regression <- calculate_gc_bias_regression(bin_gc_summary)
  
  return(gc_bias_regression)
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
  cols_to_keep <- c("sample", relative_bias_name, at_bias_name, gc_bias_name)
  gc_bias_classification <- gc_bias_loess[, ..cols_to_keep]

  # Convert bias score into fold change
  bias_cols <- grep("^bias_", names(gc_bias_classification), value = TRUE)
  for (col in bias_cols) {
    gc_bias_classification[, paste0(col, "_fold_change") := 2^get(col) - 1]
  }

  gc_bias_classification[, paste0(relative_bias_name, "_failure_threshold_at") := relative_bias_cutoff_failure_lower]
  gc_bias_classification[, paste0(relative_bias_name, "_failure_threshold_gc") := relative_bias_cutoff_failure_upper]
  gc_bias_classification[, paste0(relative_bias_name, "_warning_threshold_at") := relative_bias_cutoff_warning_lower]
  gc_bias_classification[, paste0(relative_bias_name, "_warning_threshold_gc") := relative_bias_cutoff_warning_upper]
  gc_bias_classification[, paste0(at_bias_name, "_failure_threshold") := cutoff_failure_at]
  gc_bias_classification[, paste0(gc_bias_name, "_failure_threshold") := cutoff_failure_gc]

  gc_bias_classification[, bias_type := fcase(
    get(relative_bias_name) >= relative_bias_cutoff_failure_upper, "GC biased",
    get(relative_bias_name) <= relative_bias_cutoff_failure_lower, "AT biased",
    get(gc_bias_name) >= cutoff_failure_gc, "GC biased",
    get(at_bias_name) >= cutoff_failure_at, "AT biased",
    get(relative_bias_name) >= relative_bias_cutoff_warning_upper, "GC bias warning",
    get(relative_bias_name) <= relative_bias_cutoff_warning_lower, "AT bias warning",
    default = "no bias"
  )]

  cols_to_keep <- c(
    "sample",
    grep(paste0("^", relative_bias_name), names(gc_bias_classification), value = TRUE),
    grep(paste0("^", at_bias_name), names(gc_bias_classification), value = TRUE),
    grep(paste0("^", gc_bias_name), names(gc_bias_classification), value = TRUE),
    "bias_type"
  )
  gc_bias_classification <- gc_bias_classification[, ..cols_to_keep]

  return(gc_bias_classification)
}

# TODO: think about a name change to justify including
# this function in this module.
calculate_gc_bias_loess <- function(
    gc_bias_regression,
    relative_bias_name,
    at_bias_name,
    gc_bias_name) {
  # Filter and select required columns.
  gc_bias_loess <- gc_bias_regression[
    # Filter out rows where gc_percentile is not close to GC_ANCHOR or AT_ANCHOR.
    (abs(gc_percentile - GC_ANCHOR / 100) < 1e-10) | (abs(gc_percentile - AT_ANCHOR / 100) < 1e-10),
    .(sample, gc_percentile, loess_depth)
  ]

  # Multiply gc_percentile by 100.
  gc_bias_loess[, gc_percentile := gc_percentile * 100]

  # Reshape wide using dcast.
  gc_bias_loess <- dcast(
    gc_bias_loess,
    sample ~ gc_percentile,
    value.var = "loess_depth"
  )

  # Rename columns with prefix.
  setnames(
    gc_bias_loess,
    old = as.character(c(GC_ANCHOR, AT_ANCHOR)),
    new = paste0("bias_", c(GC_ANCHOR, AT_ANCHOR))
  )

  # Compute GC-to-AT relative bias score.
  gc_bias_loess[, (relative_bias_name) := 
    log2((2^get(gc_bias_name) - 1) / (2^get(at_bias_name) - 1) + 1)
  ]

  return(gc_bias_loess)
}

plot_gc_profiles <- function(gc_bias_regression, gc_bias_classification, sample_labels) {
  label <- NULL
  gc_bias_regression_w_labels <- gc_bias_regression
  # Check if sample_labels file is provided by the user
  if (!is.null(sample_labels)) {
    # Second column name is the label
    label <- colnames(sample_labels)[2]
    gc_bias_regression_w_labels <- merge(
      gc_bias_regression_w_labels,
      sample_labels,
      by = "sample",
      all.x = TRUE
    )
  }

  # Calculate y-axis limits with more aesthetic rounding
  if (identical(Y_LIM, "auto")) {
    min_loess_depth <- min(gc_bias_regression$loess_depth)
    y_min <- ifelse(min_loess_depth < 0,
                    floor(min_loess_depth * 2) / 2,
                    0)
    y_max <- ceiling(max(gc_bias_regression$loess_depth) * 2) / 2
  } else {
    y_min <- Y_LIM[1]
    y_max <- Y_LIM[2]
  }
  
  # Generate sample GC profiles plot.
  base_plot <- merge(
    gc_bias_regression_w_labels,
    gc_bias_classification[, .(sample, bias_type)],
    by = "sample",
    all.x = TRUE
  ) %>%
    ggplot(aes(x = gc_percentile, y = loess_depth, color = bias_type))

  # Conditionally add the geom_line layer
  if (!is.null(label)) {
    geom_line_layer <- geom_line(aes(group = sample, linetype = .data[[label]]))
  } else {
    geom_line_layer <- geom_line(aes(group = sample))
  }

  # Add the rest of the layers
  p <- base_plot +
    geom_line_layer +
    scale_y_continuous("LOESS Depth Per GC Percentile",
                      limits = c(y_min, y_max)) +
    ggtitle("GC Content vs. Coverage by Sample") +
    scale_x_continuous("GC Content", 
                      limits = c(0, 1),
                      breaks = seq(0, 1.0, 0.2)) +
    scale_color_manual(name = "Bias Type", values = BIAS_COLORS) +
    theme_bw(base_size = 20) +
    theme(legend.position = "right")

  return(p)
}

plot_per_base_coverage <- function(all_libraries_raw_coverage) {
  # Calculate max depth to determine appropriate breaks
  max_depth <- max(all_libraries_raw_coverage$depth)
  max_power <- ceiling(log10(max_depth))

  # Generate major breaks (powers of 10)
  major_breaks <- 10^seq(0, max_power)

  # Generate minor breaks for each decade
  minor_breaks <- unlist(lapply(seq(0, max_power-1), function(power) {
    seq(10^power, 10^(power+1), length.out = 10)
  }))

  # Truncate sample names if they're too long
  all_libraries_raw_coverage$sample_short <- ifelse(
    nchar(all_libraries_raw_coverage$sample) > 15,
    paste0(substr(all_libraries_raw_coverage$sample, 1, 12), "..."),
    all_libraries_raw_coverage$sample
  )

  p <- ggplot(all_libraries_raw_coverage, aes(sample_short, depth)) +
    geom_boxplot(varwidth = TRUE) +
    scale_y_continuous("Coverage", 
                      trans = "log10", 
                      limits = c(1, NA),
                      breaks = major_breaks,
                      minor_breaks = minor_breaks) +
    scale_x_discrete("Sample") +
    ggtitle("Per-Base Coverage") +
    theme_bw(base_size = 20) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5))
  return(p)
}

main <- function(
    bam_coverage_directory,
    reference_gc_content_file,
    sample_labels_csv,
    outdir) {
  # Create output directory if it doesn't exist
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  # Set variable names.
  relative_bias_name <- str_glue("bias_{AT_ANCHOR}vs{GC_ANCHOR}")
  at_bias_name <- str_glue("bias_{AT_ANCHOR}")
  gc_bias_name <- str_glue("bias_{GC_ANCHOR}")

  ## ------ Read files.
  ## Read GC content
  reference_gc_content <- read_gc_content_file(reference_gc_content_file)

  ## Read sample types if provided
  sample_labels <- if (sample_labels_csv != "NO_FILE") {
    read_sample_labels_csv(sample_labels_csv)
  } else {
    NULL
  }

  ## Read coverage data for each sample.
  raw_coverage_data_tables <- read_coverage_files_and_create_data_tables(
    bam_coverage_directory,
    "_intersected_coverage\\.bed$"
  )
  all_libraries_raw_coverage <- rbindlist(
    Map(function(dt, name) {
      dt[, sample := sub("_intersected_coverage$", "", name)]
      dt
    }, raw_coverage_data_tables, names(raw_coverage_data_tables))
  )

  if (DRAW_PER_BASE_COVERAGE) {
    p_per_base_coverage <- plot_per_base_coverage(all_libraries_raw_coverage)
    n_samples <- length(unique(all_libraries_raw_coverage$sample))
    ggsave(
      file.path(outdir, "per_base_coverage.png"),
      plot = p_per_base_coverage,
      height = 7.5,
      width = 0.5 * n_samples + 3,
      units = "in"
    )
  }

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
      "R script for calculating and visualizing GC biases. Required inputs are bam files",
      "in a directory and a GC content file calculated by bin/calculate_gc_content.sh.",
      "Optional inputs are a sample labels csv file and an output directory.",
      "Usage = panelGC_main.R --bam_coverage_directory <bam_coverage_directory>",
      "--reference_gc_content_file <reference_gc_content_file>",
      "--sample_labels_csv <sample_labels_csv> (optional)",
      "--outdir <outdir>"
    )
  )
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
    default = "NO_FILE"
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
    "--y_lim",
    help = "y-axis minimum and maximum.",
    default = "auto"
  )
  parser <- add_argument(
    parser,
    "--draw_per_base_coverage",
    help = "Draw per-base coverage plot."
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
  Y_LIM <<- pluck(args, "y_lim")
  if (Y_LIM != "auto") {
    Y_LIM <<- as.numeric(str_split(Y_LIM, ",")[[1]])
    if (length(Y_LIM) != 2 || !all(is.numeric(Y_LIM)) || Y_LIM[1] >= Y_LIM[2]) {
      stop("y_lim must be a comma-separated string of two numbers where the first number is less than the second.")
    }
  }
  DRAW_TREND <<- as.logical(pluck(args, "draw_trend"))
  SHOW_SAMPLES <<- as.logical(pluck(args, "show_sample_names"))
  DRAW_PER_BASE_COVERAGE <<- as.logical(pluck(args, "draw_per_base_coverage"))
  BIN_FOLDER <<- find_here()
  main(
    bam_coverage_directory,
    reference_gc_content_file,
    sample_labels_csv,
    outdir
  )
}
