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
# Minimum number of fresh libraries required for gating.
N.LIBRARY.CUTOFF <- 5

BIAS.COLORS <- c("GC biased" = "blue",
                 "AT biased" = "red",
                 "GC bias warning" = "cornflowerblue",
                 "AT bias warning" = "coral",
                 "No bias" = "darkgrey")
