plot_gc_profiles <- function(gc_bias_regression, gc_bias_classification) {
 # Maximum of y-axis.
 y_max <- ceiling(max(pull(gc_bias_regression, normalized_depth)))
 # Generate library GC profiles plot.
 p <- gc_bias_regression %>%
  left_join(select(gc_bias_classification,
                   "library", "Bias_Type"),
              by = "library", multiple = TRUE) %>%
  # mutate(FFPE = ifelse(str_detect(`Nature Of Analyte`, "FFPE"), TRUE, FALSE),
  #        FFPE = ifelse(is.na(FFPE), FALSE, FFPE)) %>%
  ggplot(aes(x = gc_bin, y = normalized_depth, color = Bias_Type)) +
  geom_smooth(aes(group = library),
              method = "loess", formula = y ~ x, se = FALSE, fullrange = TRUE) +
  scale_y_continuous("Relative Depth Per GC Bin") +
  ggtitle(str_glue(
   "Probe GC Content vs. Coverage by Library")) +
  coord_cartesian(xlim = c(0.0, 1.0), ylim = c(0, y_max)) +
  scale_x_continuous("Probe GC", breaks = seq(0, 1.0, 0.2)) +
  scale_color_manual(values = BIAS.COLORS) +
  theme_bw(base_size = 20) +
  theme(legend.position = "right")
 return(p)
}

