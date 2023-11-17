#! /usr/bin/env Rscript

source("new_R_modules/01_load_libraries_and_constants.R")
source("new_R_modules/02_load_data.R")
source("new_R_modules/03_data_manipulators.R")
source("new_R_modules/04_gc_bias_regression.R")
source("new_R_modules/05_gc_bias_classification.R")
source("new_R_modules/06_gc_bias_plot.R")

# Set variable names.
relative_bias_name <- str_glue("Bias_{GC.LOWER.ANCHOR}vs{GC.UPPER.ANCHOR}")
at_bias_name <- str_glue("Bias_{GC.LOWER.ANCHOR}")
gc_bias_name <- str_glue("Bias_{GC.UPPER.ANCHOR}")

## development hardcoded arguments

probe_bed_file <- "/projects/clingenetics/mgoktas_dev/other-tix/GC_PANEL_PROJECT/panelGC/test_files/CCG.oncohcp.dna.probes.GRCh37.idt.2020.1.bed"
bam_coverage_directory <- "/projects/clingenetics/mgoktas_dev/other-tix/GC_PANEL_PROJECT/panelGC/work/bc/d66b111ac7446613bde4201d56365d/"
reference_gc_content_file <- "work/fd/bd48da3462788fd9b6c157a92191ab/extracted_sequences_GC_content.txt"
outdir <- "/projects/clingenetics/mgoktas_dev/other-tix/GC_PANEL_PROJECT/panelGC/test_files/test_outdir/"

## ------ Read files.
## Read probes.
probes <- read_bed_file(probe_bed_file)

## Read GC content
reference_gc_content <- read_gc_content_file(reference_gc_content_file)
## Add probe name/identifier to reference_gc_content.
reference_gc_content <- reference_gc_content %>% left_join(
	probes,
	by = c("chromosome" = "chromosome",
				 "start_pos" = "probe_start_pos",
				 "end_pos" = "probe_end_pos")
)

## Read coverage data for each library.
raw_coverage_tibbles <- read_coverage_files_and_create_tibbles(
	bam_coverage_directory, 
	"_intersected_coverage\\.bed$"
	)
all_libraries_raw_coverage <- raw_coverage_tibbles %>%
	imap(~ mutate(.x, library = .y)) %>%
	bind_rows()
## Add probe name/identifier to reference_gc_content.
all_libraries_raw_coverage <- all_libraries_raw_coverage %>% left_join(
		probes,
		by = c("chromosome" = "chromosome",
					 "probe_start_pos" = "probe_start_pos",
					 "probe_end_pos" = "probe_end_pos")
	)

all_libraries_raw_coverage <- add_offset_to_start_position_and_drop_columns(
	all_libraries_raw_coverage
	)

gc_bias_regression_table <- get_gc_bias_regression_table(
	all_libraries_raw_coverage,
	reference_gc_content)

# Write result to "gc_bias_loess_regression.tsv"

cutoffs_failure <- compute_cutoffs(RELATIVE.FOLD.CHANGE.FAILURE.CUTOFF)
cutoff_failure_upper <- cutoffs_failure$upper
cutoff_failure_lower <- cutoffs_failure$lower

cutoffs_warning <- compute_cutoffs(RELATIVE.FOLD.CHANGE.WARNING.CUTOFF)
cutoff_warning_upper <- cutoffs_warning$upper
cutoff_warning_lower <- cutoffs_warning$lower

gc_bias_classification_table <- get_gc_bias_classification_table(
	gc_bias_regression_table,
	cutoff_failure_upper, cutoff_failure_lower,
	cutoff_warning_upper, cutoff_warning_lower)
print(gc_bias_classification_table)


p <- plot_gc_profiles(gc_bias_regression_table, gc_bias_classification_table)
ggsave("my_plot.png", plot = p, height = 7.5, width = 10, units = "in")


calculate_pass_summary <- function(data, status_column, pass_value){
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

# GC bia pass ratios.
pass_summary_gc_bias <- calculate_pass_summary(
	gc_bias_classification_table,
	"Bias_Type", 
	"No bias"
	)
print("{n_pass} {n_total} {pass_ratio}")
print(pass_summary_gc_bias)

# Panel GC bias summary.
# TODO: Create and print summary.





