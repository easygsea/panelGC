#!/bin/bash

input_extracted_sequences=$1
output_gc_content=$2
paste -sd '\t\n' $input_extracted_sequences \
  | sed 's/^>//' \
  | awk -F'\t' 'BEGIN {OFS="\t"; print "region\tGC"} \
  {total_count=length($2); gsub(/A|T|a|t/, "", $2); gc_count=length($2); \
  gc_content=sprintf("%.15g",gc_count/total_count); print $1, gc_content}' \
  > $output_gc_content
