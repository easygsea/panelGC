#!/bin/bash

window_size=$1
input_extracted_sequences=$2
output_gc_content=$3

awk -F'\t' -v window_size=$window_size 'BEGIN {OFS="\t"; print "region\tregion_name\twindow_start\twindow_end\tGC"} \
  {
    region_name = $1
    region = $2
    seq = toupper($3)
    seq_len = length(seq)

    if (window_size == 0) {
      # Calculate GC content for entire sequence when window_size is 0
      gsub(/A|T/, "", seq)
      gc_count = length(seq)
      gc_content = sprintf("%.15g", gc_count/seq_len)
      print region, region_name, 1, seq_len, gc_content
      next
    }

    if (seq_len < window_size) {
      # Skip sequences shorter than window_size
      next
    }
    
    # For each position in the sequence
    for (i = 0; i <= seq_len - window_size; i++) {
      # Calculate window start and end positions
      start = i
      end = start + window_size
      
      # Extract window sequence
      window = substr(seq, start + 1, window_size)
      
      # Calculate GC content for this window
      gsub(/A|T/, "", window)
      gc_count = length(window)
      gc_content = sprintf("%.15g", gc_count/window_size)
      
      # Print result
      print region, region_name, start + 1, end, gc_content
    }
  }' \
  $input_extracted_sequences > $output_gc_content