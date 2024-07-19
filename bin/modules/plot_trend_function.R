BIAS_COLORS <- c(
  "GC biased" = "blue",
  "AT biased" = "red",
  "GC bias warning" = "cornflowerblue",
  "AT bias warning" = "coral",
  "no bias" = "darkgrey"
)

extract_cutoff <- function(
    gc_bias_classification,
    column_name,
    suffix,
    type = "failure"
  ) {
  gc_bias_classification %>%
    pull(str_glue("{column_name}_{type}_threshold{suffix}")) %>%
    sort() %>%
    first()
}

plot_geom_vline <- function(gc_bias_classification, bias_name, bias_cutoff) {
  geom_vline(
    data = filter(gc_bias_classification, Bias_Category == bias_name),
    aes(xintercept=bias_cutoff),
    linetype = "dashed", linewidth = 0.2, color = "black"
  )
}

plot_gc_trend <- function(gc_bias_classification, show_sample_names) {
  bias_names <- colnames(gc_bias_classification)
  bias_names <- sort(bias_names[str_detect(bias_names, "^bias_.*\\d$")])

  at_bias_cutoff <- extract_cutoff(gc_bias_classification, bias_names[1], "")
  gc_bias_cutoff <- extract_cutoff(gc_bias_classification, bias_names[3], "")
  relative_bias_failure_cutoff_at <- extract_cutoff(
    gc_bias_classification, bias_names[2], "_at")
  relative_bias_failure_cutoff_gc <- extract_cutoff(
    gc_bias_classification, bias_names[2], "_gc")
  relative_bias_warning_cutoff_at <- extract_cutoff(
    gc_bias_classification, bias_names[2], "_at", type = "warning")
  relative_bias_warning_cutoff_gc <- extract_cutoff(
    gc_bias_classification, bias_names[2], "_gc", type = "warning")
  
  gc_bias_classification <- gc_bias_classification %>%
    pivot_longer(
      cols = bias_names,
      names_to = "Bias_Category",
      values_to = "Bias_Score"
    )
  
  p <- gc_bias_classification %>%
    ggplot(aes(x = Bias_Score, y = sample, color = bias_type)) +
    geom_point() +
    plot_geom_vline(gc_bias_classification, bias_names[1], at_bias_cutoff) +
    plot_geom_vline(gc_bias_classification, bias_names[3], gc_bias_cutoff) +
    plot_geom_vline(gc_bias_classification, bias_names[2], relative_bias_failure_cutoff_at) +
    plot_geom_vline(gc_bias_classification, bias_names[2], relative_bias_failure_cutoff_gc) +
    plot_geom_vline(gc_bias_classification, bias_names[2], relative_bias_warning_cutoff_at) +
    plot_geom_vline(gc_bias_classification, bias_names[2], relative_bias_warning_cutoff_gc) +
    scale_color_manual(values = BIAS_COLORS) +
    labs(x = "Bias Score", y = "", color = "Bias Type") +
    theme_bw() +
    facet_grid(~ Bias_Category, scales = "free_y", space = "free")
  
  if(show_sample_names){
    p <- p +
      theme(
        strip.text.y.right = element_text(angle = 0)
      )
  }else{
    p <- p +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }
  
  return(p)
}