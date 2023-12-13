calculate_gc_bias_regression <- function(bin_gc_summary) {
  gc_bias_regression <- bin_gc_summary %>%
    group_by(library) %>%
    arrange(library, gc_bin) %>%
    nest() %>%
    mutate(Loess = purrr::map(data, function(x) {
      stats::loess(normalized_depth ~ gc_bin, span = 0.75, data = x) %>%
        stats::predict(gc_bin = unique(gc_bin))
    })) %>%
    unnest(cols = c(data, Loess))
  return(gc_bias_regression)
}

get_gc_bias_regression_table <- function(
    raw_bam_readcount_intersected_probes,
    gc_content) {
  mean_library_depths <- raw_bam_readcount_intersected_probes %>%
    select(library, chromosome, position, depth) %>%
    distinct() %>%
    group_by(library) %>%
    summarise(mean_depth = mean(depth))

  # Compute the median depth of each probe per library.
  gc_summary <- raw_bam_readcount_intersected_probes %>%
    inner_join(gc_content,
      by = c("probe_name" = "probe_name"), ## TODO: remove 'by'
      multiple = "all"
    ) %>%
    group_by(probe_name, GC, library) %>%
    summarise(depth_median = median(depth))

  # Compute the median depth across probes within a GC percentile bin.
  bin_gc_summary <- gc_summary %>%
    group_by(library) %>%
    mutate(gc_bin = floor(GC * 100) / 100) %>%
    group_by(library, gc_bin) %>%
    summarise(depth_median_median = median(depth_median))

  # Normalization by mean library depth.
  bin_gc_summary <- left_join(
    bin_gc_summary, mean_library_depths,
    by = "library"
  ) %>%
    mutate(normalized_depth = log2(
      depth_median_median / mean_depth + 1
    ))

  gc_bias_regression <- calculate_gc_bias_regression(bin_gc_summary)
}
