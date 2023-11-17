# compute_gc_bias_cutoff <- function(fold_change_cutoff, upper){
#  if(upper){
#   bias_cutoff <- log2(fold_change_cutoff + 1)
#  }else{
#   bias_cutoff <- log2(1/fold_change_cutoff + 1)
#  }
# }

compute_bias_cutoff <- function(fold_change_cutoff, is_upper){
 bias_factor <- if (is_upper) fold_change_cutoff else 1/fold_change_cutoff
 return(log2(bias_factor + 1))
}

compute_cutoffs <- function(cutoff){
 upper <- compute_bias_cutoff(cutoff, TRUE)
 lower <- compute_bias_cutoff(cutoff, FALSE)
 return(list(upper = upper, lower = lower))
}

get_gc_bias_classification_table <- function(gc_bias_regression,
                                             relative_bias_cutoff_failure_upper,
                                             relative_bias_cutoff_failure_lower,
                                             relative_bias_cutoff_warning_upper,
                                             relative_bias_cutoff_warning_lower){
 gc_bias_loess <- calculate_gc_bias_loess(gc_bias_regression)
 print(gc_bias_loess)
 gc_bias_classification <- classify_gc_bias(
  gc_bias_loess,
  relative_bias_cutoff_failure_upper, relative_bias_cutoff_failure_lower,
  relative_bias_cutoff_warning_upper, relative_bias_cutoff_warning_lower)
 return(gc_bias_classification)
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

#TODO: think about a name change to justfiy including this function in this module.
calculate_gc_bias_loess <- function(gc_bias_regression){
 gc_bias_loess <- gc_bias_regression %>%
  filter(gc_bin == GC.UPPER.ANCHOR/100 | gc_bin == GC.LOWER.ANCHOR/100) %>%
  select(library, gc_bin, Loess) %>%
  mutate(gc_bin = gc_bin * 100) %>%
  pivot_wider(names_from = gc_bin, values_from = Loess,
              names_prefix = "Bias_") %>%
  # GC-to-AT relative bias is computed as
  # log2(HQ_Depth_GCupper/HQ_Depth_GClower + 1)
  mutate(!!relative_bias_name :=
          log2((2 ** get(gc_bias_name) - 1) /
                (2 ** get(at_bias_name) - 1) + 1))
 return(gc_bias_loess)
}
