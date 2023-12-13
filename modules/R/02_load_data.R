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
  manipulated_data <- data %>%
    separate(region, into = c("chromosome", "positions"), sep = "[:]") %>%
    separate(positions, into = c("start_pos", "end_pos"), sep = "-") %>%
    mutate(
      start_pos = as.integer(start_pos),
      end_pos = as.integer(end_pos)
    )

  return(manipulated_data)
}

read_bed_file <- function(file_path) {
  message("in read bed file")
  # Replace 'your_file.bed' with the path to your BED file
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
  message("leaving read bed file")
  # Convert to tibble
  return(bed_data)
}
