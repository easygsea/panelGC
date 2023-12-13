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
# Exclusive minimum fraction of libraries showing GC bias/ONCOCNV
# calls that signals a potential issue for a panel.
# TODO: This variable is not used
PANEL_GC_BIAS_PASS_RATIO_CUTOFF <- 0.7
# Minimum number of fresh libraries required for gating.
# TODO: This variable is not used
N_LIBRARY_CUTOFF <- 5

BIAS_COLORS <- c(
  "GC biased" = "blue",
  "AT biased" = "red",
  "GC bias warning" = "cornflowerblue",
  "AT bias warning" = "coral",
  "no bias" = "darkgrey"
)
