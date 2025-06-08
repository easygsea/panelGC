#!/bin/bash

input_extracted_sequences=$1
output_gc_content=$2
awk -F'\t' 'BEGIN {OFS="\t"; print "region\tname\tGC"} \
  {total_count=length($3); gsub(/A|T|a|t/, "", $3); gc_count=length($3); \
  gc_content=sprintf("%.15g",gc_count/total_count); print $2, $1, gc_content}' \
  $input_extracted_sequences > $output_gc_content