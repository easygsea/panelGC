# Project Name - R Modules

## Overview
This directory contains separate R modules used during the development phase of the project. Each module is a standalone file that represents a distinct component or functionality of the project.

## File Naming Convention
Each R file in this directory follows a specific naming convention:
- Files are prefixed with a two-digit number (e.g., `01_`, `02_`, etc.) to denote the order of execution.
- The prefix is followed by a descriptive name, which indicates the module's purpose (e.g., `01_load_libraries_and_constants.R`).

## Development Process
During the development phase, developers are encouraged to work on separate files. This approach facilitates better maintenance and understanding of individual components of the project.

## Deployment Preparation
Before deploying the project, it is essential to concatenate these individual R modules into a single script. This ensures that the deployment package is compact and easy to manage.

### Concatenation Script
To concatenate these files, use the provided shell script `concat_R_modules_and_move_to_bin.sh`. This script automates the process of combining the R modules into a single file.

#### Script Details
The script performs the following actions:
1. **Find Modules**: Locates all R files in the current directory that match the pattern `./[0-9][0-9]_*.R`.
2. **Sort and Concatenate**: Sorts these files in ascending order and concatenates them into a single file.
3. **Remove Unnecessary Lines**: Filters out lines containing `source` and `Rscript` commands, as these are typically not needed in the combined script.
4. **Output File**: The final concatenated script is saved as `panelGC_main.R` in the `../../bin/` directory.

#### Usage
Run the following command in the terminal within this directory:
```bash
./concat_R_modules_and_move_to_bin.sh
