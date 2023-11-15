#! /usr/bin/env Rscript
# generate_gc_biases.R
# Created by Jean Cheng, Jill Slind, Jelena Sihvonen, David Cohn on 2023-08-21.
# Copyright (c) 2023 Canada's Michael Smith Genome Sciences Centre.
# All rights reserved.
# The script is applicable for non-Prostate panels only.

# Before executing the program, please run the following command to configure
# the appropriate R environment:
# source <STANDALONE_UTILITIES_ROOT>/R_utilities/R_profile.sh

libpath <- "/projects/clingenetics/software/tools/R/R-3.6.1/library"
.libPaths(libpath)

message("Compute and classify GC biases. Loading packages and functions ...")

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
GC.UPPER.ANCHOR <- 75
GC.LOWER.ANCHOR <- 25
RELATIVE.FOLD.CHANGE.FAILURE.CUTOFF <- 2
RELATIVE.FOLD.CHANGE.WARNING.CUTOFF <- 1.5
# Exlusive minimum fraction of libraries showing GC bias/ONCOCNV calls that signals a
# potential issue for a panel.
PANEL.GC.BIAS.PASS.RATIO.CUTOFF <- 0.7
PANEL.ONCOCNV.PASS.RATIO.CUTOFF <- 0.5
# Minimum number of fresh libraries required for gating.
N.LIBRARY.CUTOFF <- 5
# Minimum number of copy calls that signals a potential issue for a library.
COPY.CHANGE.NUMBER.FAIL.CUTOFF <- 10
COPY.CHANGE.NUMBER.WARNING.CUTOFF <- 5
# Minimum fraction of copy calls that signals a potential issue for a library.
COPY.CHANGE.RATIO.FAIL.CUTOFF <- 0.1
COPY.CHANGE.RATIO.WARNING.CUTOFF <- 0.05
# Panel names in analysis folders.
HCP.PANEL.NAME.IN.ANALYSIS.FOLDER <- "hereditary"
MYELOID.PANEL.NAME.IN.ANALYSIS.FOLDER <- "Myeloid"
# Panels used for pool-level GC bias evaluation.
INDICATION.PANELS <- c(HCP.PANEL.NAME.IN.ANALYSIS.FOLDER,
                       MYELOID.PANEL.NAME.IN.ANALYSIS.FOLDER)
## Regular panel deployment directory to be used after probes GC files merged
## into corresponding master branches.
# PANEL.DEPLOYMENT.DIRECTORY <- file.path("", "projects", "CDG_production",
#                                         "deployments", "panels")
## Temporary panel deployment directory for testing before probes GC files merged.
PANEL.DEPLOYMENT.DIRECTORY <- file.path(
  "", "projects", "clingenetics", "CCG_Dev", "CCG_Dev_jcheng", "clingen",
  "CLINGEN-8850.library-wise_gc_bias_gate", "CLINGEN-8851.place_probe_gc_txt")
# Production analysis directory, ensuring that "--write_to_panel_directory" is
# used when "--outdir" is configured to refer to a production pool directory.
PRODUCTION.ANALYSIS.DIRECTORY <- file.path(
  "", "projects", "CDG_production", "analysis")
BIAS.COLORS <- c("GC biased" = "blue",
                 "AT biased" = "red",
                 "GC bias warning" = "cornflowerblue",
                 "AT bias warning" = "coral",
                 "No bias" = "darkgrey")

relative_bias_name <- str_glue("Bias_{GC.LOWER.ANCHOR}vs{GC.UPPER.ANCHOR}")
at_bias_name <- str_glue("Bias_{GC.LOWER.ANCHOR}")
gc_bias_name <- str_glue("Bias_{GC.UPPER.ANCHOR}")

# Functions ----
# * Functions: Helper ----
annotate_fail_color <- function(value){
  ifelse(isFALSE(value) | value == "FALSE" | value == "Fail",
         sprintf("{color:#DE350B}%s{color}", value), value)
}
load_or_write_table <- function(filename, overwrite, analysis_function){
  filename_exists <- file.exists(filename)
  if(filename_exists & !overwrite){
    data <- read_tsv(filename, col_names = TRUE, col_types = cols())
    message(str_glue(
      "{filename} already exists. Skipped writing. Data loaded from file."))
  }else{
    data <- analysis_function
    write.table(
      data, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
    if(filename_exists){
      message(str_glue("Overwritten to {filename}"))
    }else{
      message(str_glue("Written to {filename}"))
    }
  }
  return(data)
}

print_jira_table <- function(data) {
  # Convert everything to characters, so e.g. dates and factors print properly.
  data.as.character <- data %>%
    mutate(across(everything(), as.character))
  c(
    paste('||', paste0(colnames(data.as.character), collapse = ' || '), '||'),
    sapply(
      1:nrow(data.as.character),
      function(row) {
        return(
          paste('|', paste0(data.as.character[row,], collapse = ' | '), '|'))}),
    '') %>%
    paste(collapse = '\n') %>%
    cat()
}
# * Functions: Data loading ----
read_submission_worksheet <- function(file_name) {
  raw_sheet <- suppressMessages(
      read_excel(file_name, col_names = FALSE))

  # Extract the batch from the worksheet header
  row_containing_batch_id <- pluck(raw_sheet, "...1") %>% str_which("Batch ID")
  batch <- raw_sheet[[row_containing_batch_id, 2]]

  # Read only the rectangular portion of the submission worksheet as a tibble
  col_containing_sample_id <- grep("Sample ID", raw_sheet)
  header_row <- raw_sheet %>%
    pull(col_containing_sample_id) %>%
    str_which("Sample ID")

  submission_worksheet <- raw_sheet %>%
    slice(header_row:n()) %>%
    set_names(as.character(slice(., 1))) %>%
    slice(-1) %>%
    mutate("Batch" = batch)

  return(submission_worksheet)
}

read_in_bam_readcount_all_sites_file <- function(bam_readcount_all_sites_file,
                                                 probes_bed_file,
                                                 negative_controls){
  # Read in depth per position per library.
  raw_bam_readcount <- read_tsv(bam_readcount_all_sites_file, col_types = cols(
    "Pool" = col_character(), "Library" = col_character(),
    "Target" = col_character(), "Chrom" = col_character(),
    "Position" = col_integer(), "HQ_Depth" = col_integer())) %>%
    filter(!Library %in% negative_controls)
  # Load probes genomic coordinates in BED format.
  probes_bed <- rtracklayer::import(probes_bed_file)
  # Intersect per-position depth to probes BED.
  raw_bam_readcount_intersected_probes <- raw_bam_readcount %>%
    makeGRangesFromDataFrame(
      seqnames.field = "Chrom", start.field = "Position",
      end.field = "Position", keep.extra.columns = TRUE) %>%
    mergeByOverlaps(probes_bed) %>%
    as.data.table() %>%
    select(
      "Chrom" = "..seqnames",
      "Position" = "..start",
      "HQ_Depth",
      "Pool",
      "Library",
      "Probe" = "name")
  return(raw_bam_readcount_intersected_probes)
}

# * Functions: GC bias analysis ----
calculate_mean_library_depth <- function(raw_bam_readcount_intersected_probes){
  # Calculate the mean depth across all targeted and probes bases per library.
  mean_library_depths <- raw_bam_readcount_intersected_probes %>%
    select(Pool, Library, Chrom, Position, HQ_Depth) %>%
    distinct() %>%
    group_by(Pool, Library) %>%
    summarise(Mean_HQ_Depth = mean(HQ_Depth))
  return(mean_library_depths)
}

calculate_gc_bias_regression <- function(bin_gc_summary){
  gc_bias_regression <- bin_gc_summary %>%
    group_by(Pool, Library) %>%
    arrange(Library, GC_Bin) %>%
    nest() %>%
    mutate(Loess = purrr::map(data, function(x)
      stats::loess(Normalized_HQ_Depth ~ GC_Bin, span = 0.75, data = x) %>%
        stats::predict(GC_Bin = unique(GC_Bin)))) %>%
    unnest(cols = c(data, Loess))
  return(gc_bias_regression)
}

calculate_gc_bias_loess <- function(gc_bias_regression){
  gc_bias_loess <- gc_bias_regression %>%
    filter(GC_Bin == GC.UPPER.ANCHOR/100 | GC_Bin == GC.LOWER.ANCHOR/100) %>%
    select(Pool, Library, GC_Bin, Loess) %>%
    mutate(GC_Bin = GC_Bin * 100) %>%
    pivot_wider(names_from = GC_Bin, values_from = Loess,
                names_prefix = "Bias_") %>%
    # GC-to-AT relative bias is computed as
    # log2(HQ_Depth_GCupper/HQ_Depth_GClower + 1)
    mutate(!!relative_bias_name :=
             log2((2 ** get(gc_bias_name) - 1) /
                    (2 ** get(at_bias_name) - 1) + 1))
  return(gc_bias_loess)
}

classify_gc_bias <- function(gc_bias_loess,
                             relative_bias_cutoff_failure_upper,
                             relative_bias_cutoff_failure_lower,
                             relative_bias_cutoff_warning_upper,
                             relative_bias_cutoff_warning_lower) {
  gc_bias_classification <- gc_bias_loess %>% pivot_longer(
    cols = matches("^Bias_.*[0-9]$"),
    names_to = "Bias_Category",
    values_to = "Bias_Score") %>%
    filter(Bias_Category == relative_bias_name) %>%
    mutate(
      Bias_Type = case_when(
        Bias_Score >= relative_bias_cutoff_failure_upper ~ "GC biased",
        Bias_Score <= relative_bias_cutoff_failure_lower ~ "AT biased",
        Bias_Score >= relative_bias_cutoff_warning_upper ~ "GC bias warning",
        Bias_Score <= relative_bias_cutoff_warning_lower ~ "AT bias warning",
        TRUE ~ "No bias"
      )) %>%
    # Convert Bias_Score into fold change.
    mutate(
      Bias_Score_Fold_Change = 2 ** Bias_Score - 1, .after = "Bias_Score")
  return(gc_bias_classification)
}

# * Functions: ONCOCNV statistics ----
extract_oncocnv_files <- function(panel_directory){
  # List all .summary.txt files.
  oncocnv_files <- system(str_glue(
    "find {panel_directory} -maxdepth 2 -mindepth 2 -name '*.summary.txt'"),
    intern = TRUE)
  return(oncocnv_files)
}

summarise_library_oncocnv_statistics <- function(oncocnv_file){
  # Load ONCOCNV summary and summarize non-neutral autosomal calls.
  oncocnv_statistics <- read_tsv(oncocnv_file, col_types = cols_only(
    "chr" = col_character(), "copy.number" = col_number())) %>%
    filter(str_detect(chr, "[0-9]$")) %>%
    mutate(copy.called = ifelse(copy.number != 2, TRUE, FALSE)) %>%
    pull(copy.called)
  # Check if all TRUE or all FALSE.
  if(all(oncocnv_statistics)){
    oncocnv_statistics <- c(length(oncocnv_statistics), 0)
    names(oncocnv_statistics) <- c("TRUE", "FALSE")
  }else if(all(!oncocnv_statistics)){
    oncocnv_statistics <- c(0, length(oncocnv_statistics))
    names(oncocnv_statistics) <- c("TRUE", "FALSE")
  }else{
    oncocnv_statistics <- table(oncocnv_statistics)
  }
  oncocnv_statistics <- tibble(Library = basename(dirname(oncocnv_file)),
                               CopyChange = names(oncocnv_statistics),
                               Count = as.numeric(oncocnv_statistics)) %>%
    pivot_wider(names_from = CopyChange, values_from = Count,
                names_prefix = "CopyChange") %>%
    relocate(CopyChangeTRUE, .after = "Library") %>%
    mutate(Ratio_CopyChangeTRUE = CopyChangeTRUE /
             (CopyChangeTRUE + CopyChangeFALSE))
  return(oncocnv_statistics)
}

summarise_panel_oncocnv_statistics <- function(oncocnv_statistics){
  oncocnv_statistics_summary <- oncocnv_statistics %>%
    rowwise() %>%
    mutate(Status = case_when(
      # Add a failure note if exceeding the failing amount of copy calls.
      CopyChangeTRUE >= COPY.CHANGE.NUMBER.FAIL.CUTOFF |
        Ratio_CopyChangeTRUE >= COPY.CHANGE.RATIO.FAIL.CUTOFF ~ "Fail",
      # Add a warning note if exceeding the warning amount of copy calls.
      CopyChangeTRUE >= COPY.CHANGE.NUMBER.WARNING.CUTOFF |
        Ratio_CopyChangeTRUE >= COPY.CHANGE.RATIO.WARNING.CUTOFF ~ "Warning",
      TRUE ~ "Pass"
    ))
  return(oncocnv_statistics_summary)
}

# * Functions: Wrappers ----
create_manifest <- function(pool_directory, uuids){
  # List panel analysis directories.
  panel_directories <- list.files(
    pool_directory, pattern = "^(Onco|Mye|here)", full.names = FALSE)
  # The analyzed panels in the pool.
  panels <- unique(sub("_.+$", "", panel_directories))
  # Extract the latest analysis for each panel.
  panel_directories <- sapply(panels, function(panel){
    panel_analyses <- sort(panel_directories[str_starts(panel_directories, panel)])
    # Extract the latest panel analysis.
    latest_panel_analysis <- last(panel_analyses)
    # Use the specified UUID if supplied.
    if(!is.na(uuids)){
      # Analyses UUIDs for the panel.
      panel_uuids <- sub("^.*_", "", panel_analyses)
      # Use the latest supplied UUID if found within the analyses UUIDs.
      uuid_index <- which(panel_uuids %in% uuids) %>% last()
      if(is.na(uuid_index)){
        message(str_glue("{panel}: {last(panel_uuids)} (latest)."))
      }else{
        # Update the panel analysis to return according to the user input.
        latest_panel_analysis <- pluck(panel_analyses, uuid_index)
        message(str_glue(
          "{panel}: {pluck(panel_uuids, uuid_index)} (user-defined)."))
      }
    }else{
      message(str_glue("{panel}: {sub('^.*_', '', latest_panel_analysis)} (latest)"))
    }
    return(latest_panel_analysis)
  })
  # Create the manifest.
  manifest <- tibble(
    pool = basename(pool_directory),
    panel = sub("_.+$", "", panel_directories),
    panel_directory = file.path(pool_directory, panel_directories),
    group_analysis_id = sub("^.*_", "", basename(panel_directory)),
    panel_index = file.path(panel_directory, str_glue("{pool}.indices")),
    panel_bam_readcount = file.path(
      panel_directory, str_glue("{pool}.bam_readcount_all_sites.csv")),
    probe_bed = list.files(panel_directory, pattern = "dna.probes.*.bed$",
                           full.names = TRUE)) %>%
    rowwise() %>%
    mutate(probe_gc = list.files(
      PANEL.DEPLOYMENT.DIRECTORY, pattern = str_glue("^{tolower(panel)}"),
      recursive = FALSE, full.names = TRUE),
      probe_gc = list.files(
        file.path(probe_gc, "resources", "probes"), pattern = ".GC.txt$",
        recursive = FALSE, full.names = TRUE))
  return(manifest)
}

extract_nature_of_analyte_from_submission <- function(file_names){
  map_dfr(file_names, read_submission_worksheet) %>%
    select("Library", "Nature Of Analyte")
}

add_nature_of_analyte_column <- function(data_source, ffpe_libraries){
  data_source %>%
    left_join(ffpe_libraries, by = "Library", multiple = TRUE) %>%
    relocate(`Nature Of Analyte`, .after = "Library")
}

compute_gc_bias_cutoff <- function(fold_change_cutoff, upper){
  if(upper){
    bias_cutoff <- log2(fold_change_cutoff + 1)
  }else{
    bias_cutoff <- log2(1/fold_change_cutoff + 1)
  }
}

get_gc_bias_regression_table <- function(panel_index_file,
                                         bam_readcount_all_sites_file,
                                         probes_bed_file,
                                         gc_content_file){
  # Load panel index data and extract negative controls.
  negative_controls <- read_tsv(panel_index_file, col_names = FALSE,
                                col_types = cols()) %>%
    filter(X6 == "negative control") %>% pull(X1)
  # Load raw bam read count data, intersect to probes BED, and exclude negative
  # controls.
  raw_bam_readcount_intersected_probes <- read_in_bam_readcount_all_sites_file(
    bam_readcount_all_sites_file, probes_bed_file, negative_controls)
  # Calculate mean depth per library.
  mean_library_depths <- calculate_mean_library_depth(
    raw_bam_readcount_intersected_probes)
  # Load probe GC content.
  gc_content <- read_tsv(gc_content_file, col_types = cols())
  # Compute the median depth of each probe per library.
  gc_summary <- raw_bam_readcount_intersected_probes %>%
    inner_join(gc_content, by = c("Probe" = "region"), multiple = "all") %>%
    group_by(Probe, Pool, GC, Library) %>%
    summarise(HQ_Depth_median = median(HQ_Depth))
  # Compute the median depth across probes within a GC percentile bin.
  bin_gc_summary <- gc_summary %>%
    group_by(Pool, Library) %>%
    mutate(GC_Bin = floor(GC * 100) / 100) %>%
    group_by(Pool, Library, GC_Bin) %>%
    summarise(HQ_Depth_median_median = median(HQ_Depth_median))
  # Normalization by mean library depth.
  bin_gc_summary <- left_join(
    bin_gc_summary, mean_library_depths, by = c("Pool", "Library")) %>%
    mutate(Normalized_HQ_Depth = log2(
      HQ_Depth_median_median / Mean_HQ_Depth + 1))
  # Perform LOESS regression on normalized relative depths.
  gc_bias_regression <- calculate_gc_bias_regression(bin_gc_summary)
  return(gc_bias_regression)
}

get_gc_bias_classification_table <- function(gc_bias_regression,
                                             relative_bias_cutoff_failure_upper,
                                             relative_bias_cutoff_failure_lower,
                                             relative_bias_cutoff_warning_upper,
                                             relative_bias_cutoff_warning_lower){
  gc_bias_loess <- calculate_gc_bias_loess(gc_bias_regression)
  gc_bias_classification <- classify_gc_bias(
    gc_bias_loess,
    relative_bias_cutoff_failure_upper, relative_bias_cutoff_failure_lower,
    relative_bias_cutoff_warning_upper, relative_bias_cutoff_warning_lower)
  return(gc_bias_classification)
}

save_figure <- function(figure_name, height = 7.5, width = 10, units = "in") {
  message(str_glue("Saving figure as {figure_name}.png"))
  ggsave(
    str_glue("{figure_name}.png"), height = height, width = width, units = units)
}

plot_gc_profiles <- function(pool, panel, gc_bias_regression, gc_bias_classification) {
  # Maximum of y-axis.
  y_max <- ceiling(max(pull(gc_bias_regression, Normalized_HQ_Depth)))
  # Generate library GC profiles plot.
  p <- gc_bias_regression %>%
    left_join(select(gc_bias_classification,
                     "Pool", "Library", "Nature Of Analyte", "Bias_Type"),
              by = c("Pool", "Library"), multiple = TRUE) %>%
    mutate(FFPE = ifelse(str_detect(`Nature Of Analyte`, "FFPE"), TRUE, FALSE),
           FFPE = ifelse(is.na(FFPE), FALSE, FFPE)) %>%
    ggplot(aes(x = GC_Bin, y = Normalized_HQ_Depth, color = Bias_Type)) +
    geom_smooth(aes(group = Library, linetype = FFPE),
                method = "loess", formula = y ~ x, se = FALSE, fullrange = TRUE) +
    scale_y_continuous("Relative Depth Per GC Bin") +
    ggtitle(str_glue(
      "{pool} {panel}\n",
      "Probe GC Content vs. Coverage by Library")) +
    coord_cartesian(xlim = c(0.0, 1.0), ylim = c(0, y_max)) +
    scale_x_continuous("Probe GC", breaks = seq(0, 1.0, 0.2)) +
    scale_color_manual(values = BIAS.COLORS) +
    theme_bw(base_size = 20) +
    theme(legend.position = "right")
  return(p)
}

generate_oncocnv_statistics <- function(panel_directory){
  oncocnv_files <- extract_oncocnv_files(panel_directory)
  if(length(oncocnv_files) > 0){
    panel_oncocnv_statistics <- map_dfr(
      oncocnv_files, summarise_library_oncocnv_statistics) %>%
      summarise_panel_oncocnv_statistics()
  }else{
    panel_oncocnv_statistics <- NULL
  }
  return(panel_oncocnv_statistics)
}

calculate_pass_ratio <- function(data, status_column, pass_value){
  if(nrow(data) == 0){
    pass_ratio <- NA
  }else{
    status_vector <- pull(data, status_column)
    n_total <- length(status_vector)
    n_pass <- sum(status_vector == pass_value)
    pass_ratio <- round(n_pass / n_total, digits = 3)
    pass_ratio <- str_glue("{n_pass} {n_total} {pass_ratio}")
  }
  return(pass_ratio)
}

generate_gc_biases <- function(pool,
                               panel,
                               panel_directory,
                               group_analysis_id,
                               panel_index,
                               panel_bam_readcount,
                               probe_bed,
                               probe_gc,
                               qa_directory,
                               ffpe_libraries,
                               overwrite){
  message(str_glue(
    "\n\nPerforming analysis on {pool} {panel} {group_analysis_id} ..."))
  # Create the qc subfolder if not existed yet.
  if(!dir.exists(qa_directory)){
    dir.create(qa_directory)
    message(str_glue("Created folder: {qa_directory}"))
  }
  # The file recording LOESS regression results.
  outfile_gc_bias_regression <- file.path(
    qa_directory, paste(pool, panel, group_analysis_id,
                        "gc_bias_loess_regression.tsv", sep = "."))
  # Check if LOESS regression table exists and overwrite disabled. If yes, load;
  # else compute the relative HQ depth for each GC bin per library, perform
  # LOESS regression, and obtain predicted relative HQ depth for each GC bin.
  gc_bias_regression <- load_or_write_table(
    outfile_gc_bias_regression, overwrite, get_gc_bias_regression_table(
      panel_index, panel_bam_readcount, probe_bed, probe_gc))
  # The upper and lower cutoffs for relative fold changes.
  relative_bias_cutoff_failure_upper <- compute_gc_bias_cutoff(
    RELATIVE.FOLD.CHANGE.FAILURE.CUTOFF, TRUE)
  relative_bias_cutoff_failure_lower <- compute_gc_bias_cutoff(
    RELATIVE.FOLD.CHANGE.FAILURE.CUTOFF, FALSE)
  relative_bias_cutoff_warning_upper <- compute_gc_bias_cutoff(
    RELATIVE.FOLD.CHANGE.WARNING.CUTOFF, TRUE)
  relative_bias_cutoff_warning_lower <- compute_gc_bias_cutoff(
    RELATIVE.FOLD.CHANGE.WARNING.CUTOFF, FALSE)
  # The file recording GC bias failures and warnings.
  outfile_gc_bias_classification <- file.path(
    qa_directory, paste(pool, panel, group_analysis_id,
                        "gc_bias_loess_classification.tsv", sep = "."))
  # Check if GC bias failures and warnings table exists and overwrite disabled.
  # If yes, load; else classify GC bias failure/warning per library based on
  # LOESS predictions.
  gc_bias_classification <- load_or_write_table(
    outfile_gc_bias_classification,
    overwrite,
    get_gc_bias_classification_table(
      gc_bias_regression,
      relative_bias_cutoff_failure_upper, relative_bias_cutoff_failure_lower,
      relative_bias_cutoff_warning_upper, relative_bias_cutoff_warning_lower) %>%
      add_nature_of_analyte_column(ffpe_libraries))
  # The file saving the GC profiles plot.
  outfile_gc_profiles_plot <- file.path(
    qa_directory, paste(pool, panel, group_analysis_id, "gc_profiles", sep = "."))
  # Check if GC profiles plot exists and overwrite disabled. If yes, skip; else,
  # generate the plot.
  if(file.exists(str_glue("{outfile_gc_profiles_plot}.png")) & !overwrite){
    message(str_glue("{outfile_gc_profiles_plot}.png already exists. ",
                     "Skipped writing."))
  }else{
    p <- plot_gc_profiles(pool, panel, gc_bias_regression, gc_bias_classification)
    save_figure(outfile_gc_profiles_plot)
  }
  # GC bia pass ratios.
  pass_ratio_gc_bias <- calculate_pass_ratio(
    filter(gc_bias_classification,
           !str_detect(`Nature Of Analyte`, "FFPE") | is.na(`Nature Of Analyte`)),
    "Bias_Type", "No bias")
  # Panel GC bias summary.
  panel_gc_bias_summary <- tibble(
    Pool = pool,
    Panel = panel,
    UUID = group_analysis_id,
    Total_library_count_minimum = N.LIBRARY.CUTOFF,
    GC_bias_pass_ratio = pass_ratio_gc_bias,
    GC_bias_pass_ratio_exclusive_minimum = PANEL.GC.BIAS.PASS.RATIO.CUTOFF,
    GC_bias_pass = NA,
    ONCOCNV_calls_pass_over_total = NA,
    ONCOCNV_calls_pass_ratio = NA,
    ONCOCNV_calls_pass_ratio_exclusive_minimum = PANEL.ONCOCNV.PASS.RATIO.CUTOFF,
    ONCOCNV_calls_pass = NA) %>%
    separate(GC_bias_pass_ratio, sep = " ", convert = TRUE,
             c("GC_bias_pass_count", "GC_bias_total_count", "GC_bias_pass_ratio")) %>%
    mutate(GC_bias_pass_over_total =
             str_glue("{GC_bias_pass_count}/{GC_bias_total_count}"),
           .before = "GC_bias_pass_ratio")
  # ONCOCNV statistics and GC bias summary for indicator panels.
  if(panel %in% INDICATION.PANELS){
    # Classify number of ONCOCNV calls as Fail, Warning, or Pass.
    panel_oncocnv_statistics <- generate_oncocnv_statistics(panel_directory)
    if(!is.null(panel_oncocnv_statistics)){
      # The file recording ONCOCNV statistics.
      outfile_oncocnv_statistics <- file.path(
        qa_directory, paste(pool, panel, group_analysis_id,
                            "oncocnv_statistics.tsv", sep="."))
      # Check if ONCOCNV statistics table exists and overwrite disabled. If yes,
      # load; else, write.
      load_or_write_table(
        outfile_oncocnv_statistics, overwrite, panel_oncocnv_statistics)
      # ONCOCNV pass ratio.
      pass_ratio_oncocnv <- calculate_pass_ratio(
        panel_oncocnv_statistics, "Status", "Pass")
      panel_gc_bias_summary <- panel_gc_bias_summary %>%
        mutate(ONCOCNV_calls_pass_ratio = pass_ratio_oncocnv) %>%
        separate(ONCOCNV_calls_pass_ratio, sep = " ", convert = TRUE,
                 c("ONCOCNV_calls_pass_count", "ONCOCNV_calls_total_count",
                   "ONCOCNV_calls_pass_ratio")) %>%
        mutate(ONCOCNV_calls_pass_over_total =
                 str_glue("{ONCOCNV_calls_pass_count}/{ONCOCNV_calls_total_count}"),
               .before = "ONCOCNV_calls_pass_ratio")
    }else{
      message(str_glue("ONCOCNV results are unavailable in {panel_directory}"))
    }
    # Determine whether the panel passed defined thresholds.
    panel_gc_bias_summary <- panel_gc_bias_summary %>%
      mutate(
        GC_bias_pass = ifelse(
          GC_bias_pass_ratio <= PANEL.GC.BIAS.PASS.RATIO.CUTOFF &
            GC_bias_total_count >= N.LIBRARY.CUTOFF, FALSE, TRUE),
        ONCOCNV_calls_pass = ifelse(
          is.na(ONCOCNV_calls_pass_ratio), NA,
          ifelse(ONCOCNV_calls_pass_ratio <= PANEL.ONCOCNV.PASS.RATIO.CUTOFF &
                   ONCOCNV_calls_total_count >= N.LIBRARY.CUTOFF, FALSE, TRUE)))
  }
  # The file saving panel-wise GC bias summary.
  outfile_panel_gc_bias_summary <- file.path(
    qa_directory,
    paste(pool, panel, group_analysis_id, "gc_bias_summary.tsv", sep = "."))
  # Check if GC bias summary table exists and overwrite disabled. If yes, load;
  # else, write.
  panel_gc_bias_summary <- load_or_write_table(
    outfile_panel_gc_bias_summary, overwrite, panel_gc_bias_summary)
  return(panel_gc_bias_summary)
}

evaluate_pool_gc_biases <- function(manifest, panel_gc_biases_summary){
  pool_name <- pull(panel_gc_biases_summary, Pool) %>% unique()
  # Evaluate if a panel is failing GC bias QC gate threshold:
  # hereditary: Failing both PANEL.GC.BIAS.PASS.RATIO.CUTOFF and PANEL.ONCOCNV.PASS.RATIO.CUTOFF
  # Myeloid: Failing PANEL.GC.BIAS.PASS.RATIO.CUTOFF
  # OncoPanel: "NA"
  panel_gc_biases_summary_info <- panel_gc_biases_summary %>%
    mutate(GC_bias_status = case_when(
      Panel == "OncoPanel" ~ "NA",
      Panel == HCP.PANEL.NAME.IN.ANALYSIS.FOLDER &
        !GC_bias_pass & !ONCOCNV_calls_pass ~ "Fail",
      Panel == MYELOID.PANEL.NAME.IN.ANALYSIS.FOLDER & !GC_bias_pass ~ "Fail",
      TRUE ~ "Pass"
    ))
  # Extract the failing status for each panel.
  extract_panel_status <- function(panel){
    status <- filter(panel_gc_biases_summary_info, Panel == panel) %>%
      pull(GC_bias_status)
  }
  # Evaluate whether the combination of panels' failing status indicates a
  # systemic GC bias issue within the pool.
  hcp_status <- extract_panel_status(HCP.PANEL.NAME.IN.ANALYSIS.FOLDER)
  myeloid_status <- extract_panel_status(MYELOID.PANEL.NAME.IN.ANALYSIS.FOLDER)
  combined_hcp_myeloid_status <- panel_gc_biases_summary_info %>%
    mutate(GC_bias_status_myeloid = myeloid_status) %>%
    filter(Panel == HCP.PANEL.NAME.IN.ANALYSIS.FOLDER) %>%
    mutate(GC_bias_status_hcp_myeloid = ifelse(
      (!GC_bias_pass | !ONCOCNV_calls_pass) & GC_bias_status_myeloid == "Fail",
      "Fail", "Pass")) %>%
    pull(GC_bias_status_hcp_myeloid)
  pool_fail <- hcp_status == "Fail" | combined_hcp_myeloid_status == "Fail"
  # Helper script to attach GC profile plots to JIRA.
  message("\n========Helper: Attach GC Profile Plots========")
  gc_profile_figures <- manifest %>%
    mutate(gc_profile_figure = str_glue("{pool}.{panel}.{group_analysis_id}.gc_profiles.png"),
           gc_profile_figure_path = file.path(qa_directory, gc_profile_figure),
           gc_profile_figure_thumbnail = str_glue("!{gc_profile_figure}|thumbnail!")) %>%
    select(panel, gc_profile_figure_path, gc_profile_figure_thumbnail)
  message("\nsource <PRODUCTION_PYTHON_VENV>")
  message(str_glue("attach_to_ticket.py <PROCESSING_TICKET> ",
                   "{paste0(pull(gc_profile_figures, gc_profile_figure_path), ",
                   "collapse = ' ')}"))
  message("\n=================End of Helper=================")
  # Begin JIRA output.
  message("\n===============BEGIN JIRA OUTPUT===============")
  if(pool_fail){
    message(str_glue("\n\nh2. {pool_name} shows systemic GC biases:\n\n"))
  }else{
    message(str_glue("\n\nh2. {pool_name} does not show systemic GC biases:\n\n"))
  }
  # Extract failing status.
  pool_gc_bias_status <- panel_gc_biases_summary %>%
    filter(Panel %in% INDICATION.PANELS) %>%
    select(Pool, Panel, UUID, ends_with("_pass")) %>%
    pivot_longer(ends_with("_pass"), names_to = "Indicator", values_to = "Pass") %>%
    filter(!(Panel == MYELOID.PANEL.NAME.IN.ANALYSIS.FOLDER & Indicator == "ONCOCNV_calls_pass")) %>%
    # Replace _ with whitespace for better JIRA display.
    mutate(Indicator = str_replace_all(Indicator, "_", " "),
           # Reword "pass" as "gate" to avoid confusion.
           Indicator = str_replace_all(Indicator, " pass$", " gate"),
           # Annotate text in red if failing the criterium.
           Pass = annotate_fail_color(Pass))
  # Print out tables w/ supporting information.
  message("h4. Pass status for indicator panels.")
  print_jira_table(pool_gc_bias_status)
  message("\nh4. Summary statistics.")
  panel_gc_biases_summary_info %>%
    # Add red warning colors to failures.
    mutate(across(ends_with("status"), annotate_fail_color)) %>%
    # Remove intermediate columns.
    select(-ends_with("_count")) %>%
    # Replace _ with whitespace in column names for JIRA display.
    rename_all(~str_replace_all(., "_", " ")) %>%
    print_jira_table()
  # Print out GC profile thumbnails.
  message("\nh4. GC profile plots.")
  gc_profile_figures %>%
    select(panel, gc_profile_figure_thumbnail) %>%
    pivot_wider(names_from = panel, values_from = gc_profile_figure_thumbnail) %>%
    print_jira_table() %>% message()
  # End JIRA output.
  message("\n===============END JIRA OUTPUT===============")
}

# * Functions: args ----
parse_args_function <- function(){
  parser <- arg_parser("R script for calculating and visualizing library GC
 biases for each panel and pool. Usage = generate_gc_biases.R --pool_directory
 <pool_directory> --submission_worksheet <submission_worksheets> --ffpe <ffpe_libraries>
 --outdir <outdir> --uuid <uuids> --write_to_panel_directory --overwrite")
  parser <- add_argument(parser, "--pool_directory", help = "Pool analysis
 directory, e.g. /projects/CDG_production/analysis/IX11882")
  parser <- add_argument(parser, "--submission_worksheet", help = "Submission
 worksheet(s), as defined in CLINGEN.0027. FFPE information will be extracted
 from the designated Excel datasheet(s). Please quote your input seperated by
 spaces if it has multiple entries. This option can be used in conjunction with
 --ffpe if additional FFPE libraries were pooled alongside with the libraries
 specified in the submission datasheet(s). If the submission datasheet is
 unavailable, use --ffpe to specify any FFPE libraries that were pooled, if applicable.")
  parser <- add_argument(parser, "--ffpe", help = "FFPE library IDs. Use this option
 when the submission datasheet is unavailable or additional FFPE libraries were
 pooled alongside. Please quote your input seperated by spaces if it has multiple
 entries.")
  parser <- add_argument(parser, "--outdir", help = "By default, outputs are saved
 in the current working directory under ./qa subfolder when the --outdir option
 is not specified. However, if the --outdir option is provided, outputs will be
 directed to the specified directory.")
  parser <- add_argument(parser, "--uuid", help = "Analysis UUIDs. When multiple
 analyses are available for a panel, this option allows you to specify analysis
 UUIDs to use, overriding the default use of the latest analysis of the panel.
 Please quote your input seperated by spaces if it has multiple entries.")
  parser <- add_argument(parser, "--write_to_panel_directory",
  help = "Specify that the output path (--outdir) is the same with the pool
 analysis directory (--pool_directory), that it follows the structure typically
 found in a production analysis folder for a pool, and that the intention is to
 store results within the respective panel folders under --pool_directory.",
  flag = TRUE)
  parser <- add_argument(parser, "--overwrite", help = "Please exercise caution
 when using this option, as it will overwrite existing results.", flag = TRUE)
  argv <- parse_args(parser)
  return(argv)
}

# Main program ----
main <- function(pool_directory,
                 submission_worksheet,
                 ffpe,
                 outdir,
                 uuid,
                 write_to_panel_directory,
                 overwrite){
  message("Checking parameters ...")
  # Check if valid pool ID entry.
  if(is.na(pool_directory)){
    message("\nPlease enter the pool analysis directory, using --pool_directory\n")
    q()
  }else if(!dir.exists(pool_directory)){
    message(str_glue("\n\n{pool_directory} does not exist. Please enter a valid ",
                     "pool analysis directory.\n\n"))
    q()
  }else if(!str_detect(basename(pool_directory), "IX[0-9]+$")){
    message(str_glue("\n\nThe basename of {pool_directory} does not look like ",
                     "a pool ID. Please enter a valid pool analysis directory.\n\n"))
    q()
  }
  # Verbose reminder about the output directory.
  outdir_common_reminder <- str_glue(
    "If this was not intended, please press Ctrl+C to stop the process and ",
    "then provide the correct output path using the --outdir option.\n\n")
  if(is.na(outdir)){
    message(str_glue(
      "\n\nResults will be written into the current working directory, in ",
      "./qa/. {outdir_common_reminder}"))
  }else{
    if(write_to_panel_directory){
      if(outdir == pool_directory){
        message(str_glue(
          "\n\nResults will be written into respective qa/ directory under each ",
          "panel directory in {pool_directory}. {outdir_common_reminder}"))
      }else{
        message(str_glue(
          "\n\nThe pool analysis directory (--pool_directory) and the output ",
          "path (--out_dir) are different: {pool_directory} vs {out_dir}. Please ",
          "do not use --write_to_panel_directory, or ensure that --out_dir ",
          "matches --pool_directory if the intention is to store results within ",
          "the respective panel folders under {pool_directory}."))
        q()
      }
    }else{
      if(dirname(outdir) == PRODUCTION.ANALYSIS.DIRECTORY){
        message(str_glue(
          "\n\nCAUTION!! If the intention is to write results into the production ",
          "folder, please use \"--write_to_panel_directory\" to ensure that the ",
          "results are appropriately placed within the respective panel ",
          "directories under {outdir}."))
        q()
      }else{
        message(str_glue(
          "\n\nResults will be written into {outdir}/qa/. {outdir_common_reminder}"))
      }
    }
  }
  # Check --uuid input.
  if(is.na(uuid)){
    uuids <- uuid
    message(paste0("The latest analysis will be extracted for each panel. Use ",
                   "--uuid to override this default behavior."))
  }else{
    uuids <- str_split(uuid, " ") %>% unlist()
  }
  manifest <- create_manifest(pool_directory, uuids) %>%
    mutate(qa_directory = ifelse(is.na(outdir), getwd(), outdir),
           qa_directory = file.path(qa_directory, "qa"))
  if(write_to_panel_directory){
    manifest <- manifest %>%
      mutate(qa_directory = file.path(panel_directory, "qa"))
  }
  # Annotate FFPE libraries.
  if(is.na(submission_worksheet) & is.na(ffpe)){
    message("Please use --submission_worksheet to specify full paths to the sample ",
            "submission Excel datasheet(s), and/or --ffpe to specify a list of ",
            "FFPE libraries, separated by spaces and enclosed in quotation marks.")
    q()
  }else{
    if(!is.na(submission_worksheet)){
      submission_worksheets <- str_split(submission_worksheet, " ") %>% unlist()
      if(all(grepl(".xlsx$", submission_worksheets))){
        # Extract FFPE libraries according to provided sample submission datasheets.
        if(all(file.exists(submission_worksheets))){
          ffpe_libraries <- extract_nature_of_analyte_from_submission(submission_worksheets)
        }else{
          # Report unfound sample submission files.
          unfound_submission_files <- submission_worksheets[!file.exists(submission_worksheets)]
          message(str_glue(
            "Unfound file(s). Please double check the path(s): ",
            paste0(unfound_submission_files, collapse = " ")))
          q()
        }
      }else{
        message(str_glue(
          "Your --submission_worksheet input contains non-Excel input. Please make ",
          "sure all files end with .xlsx: {submission_worksheets}"))
        q()
      }
    }
    if(!is.na(ffpe)){
      if(all(grepl("^(F|Z)[0-9]+$", ffpe))){
        ffpe <- str_split(ffpe, " ")[[1]]
        # Record provided FFPE library IDs.
        ffpe_libraries_list <- tibble(Library = ffpe,
                                      `Nature Of Analyte` = "FFPE DNA")
        if(exists("ffpe_libraries")){
          ffpe_libraries <- rbind(ffpe_libraries, ffpe_libraries_list) %>%
            distinct()
        }else{
          ffpe_libraries <- ffpe_libraries_list
        }
      }else{
        message(paste0("Please check your --ffpe input, which should be a list ",
                       "of FFPE libraries, separated by spaces and enclosed in ",
                       "quotation marks."))
        q()
      }
    }
  }
  # Attach FFPE libraries and overwrite option to the manifest file.
  manifest <- manifest %>%
    rowwise() %>%
    mutate(ffpe_libraries = list(ffpe_libraries),
           overwrite = overwrite)
  # Loop over each panel.
  panel_gc_biases_summary <- pmap_dfr(manifest, generate_gc_biases)
  evaluate_pool_gc_biases(manifest, panel_gc_biases_summary)
  message("\nAnalysis complete!")
}

if(!interactive()) {
  args <- parse_args_function()
  pool_directory <- pluck(args, "pool_directory")
  submission_worksheet <- pluck(args, "submission_worksheet")
  ffpe <- pluck(args, "ffpe")
  outdir <- pluck(args, "outdir")
  uuid <- pluck(args, "uuid")
  write_to_panel_directory <- pluck(args, "write_to_panel_directory")
  overwrite <- pluck(args, "overwrite")
  main(pool_directory, submission_worksheet, ffpe, outdir, uuid,
       write_to_panel_directory, overwrite)
}