add_offset_to_start_position_and_drop_columns <- function(tibble) {
  tibble %>%
    # Add offset to start position. Subtract one as offset starts with 1.
    mutate(position = probe_start_pos + offset - 1) %>%
    # Select the required columns
    select(library, probe_name, chromosome, position, depth, )
}
