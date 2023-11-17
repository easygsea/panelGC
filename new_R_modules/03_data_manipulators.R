add_offset_to_start_position_and_drop_columns <- function(tibble) {
  tibble %>%
   mutate(position = probe_start_pos + offset - 1) %>%  # Add offset to start position. Subtract one as offset starts with 1.
   select(library, probe_name, chromosome, position, depth, )  # Select the required columns
}
