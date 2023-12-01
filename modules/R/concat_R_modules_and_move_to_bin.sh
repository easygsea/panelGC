#!/bin/bash

# Create or overwrite the file with the Rscript shebang
echo '#! /usr/bin/env Rscript' > ../../bin/panelGC_main.R

# Concatenate the R modules into the file
find . -maxdepth 1 -regex './[0-9][0-9]_.*\.R' | sort | xargs cat | grep -v "source" | grep -v "Rscript" >> ../../bin/panelGC_main.R

# Make the script executable
chmod +x ../../bin/panelGC_main.R
