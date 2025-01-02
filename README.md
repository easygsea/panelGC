# panelGC

## Introduction
panelGC effectively quantifies and monitors GC bias in hybridization capture sequencing, identifying and flagging potential sample and/or procedural anomalies.

## Supported Architectures
amd64, arm64v8, ppc64le, s390x ([more info](https://github.com/docker-library/official-images#architectures-other-than-amd64))

## Dependencies
- Singularity (tested with 2.5.1)

Containerized dependencies (no installation required):
Nextflow 24.10.0
bedtools 2.30.0
r-base 4.3.2 argparser 0.7.1 BiocManager 1.30.22 GenomicRanges 1.54.1 rtracklayer 1.62.0 tidyverse 2.0.0

## Installation
- Install Singularity: \
Refer to official [Singularity installation guide](https://docs.sylabs.io/guides/latest/user-guide/quick_start.html)

## Usage
### Pull Singularity Image:
```bash
singularity pull docker://murathangoktas/panelgc:latest
```
### Basic Usage:
```bash
singularity exec \
  -B $(pwd):/workspace /path/to/panelgc_latest.sif \
  nextflow run /opt/panelGC/panelGC.nf \
  --bam_directory_path /workspace/bam_files/ \
  --bed_file_path /workspace/bed_file.bed \
  --fasta_file_path /workspace/fasta_file.fa \
  --out_dir /workspace/output_directory/
```
### Advanced Usage (Data Outside Working Directory):
```bash
singularity exec \
  -B $(pwd):/workspace,<data_path>:/data /path/to/panelgc_latest.sif \
  nextflow run /opt/panelGC/panelGC.nf \
  --bam_directory_path /data/bam_files/ \
  --bed_file_path /workspace/bed_file.bed \
  --fasta_file_path /workspace/fasta_file.fa \
  --out_dir /workspace/output_directory/
```
Replace the paths with your actual mounted data directories and file paths. See [Singularity documentation](https://docs.sylabs.io/guides/2.5/user-guide/bind_paths_and_mounts.html) for more information on binding paths and mounts. 

**Note:** `singularity exec` and `nextflow` provide additional options to help optimizeity the deployment of panelGC for your specific use. For more details, please refer to the official guides: [singularity exec](https://docs.sylabs.io/guides/latest/user-guide/cli/singularity_exec.html) and [nextflow CLI reference](https://www.nextflow.io/docs/latest/reference/cli.html#options).

### Parameters
- --bam_directory_path: Path to the directory containing alignment BAM files. Indices are preferred but not mandatory. Symlinks to the BAM and index files are valid.
- --bed_file_path: Path to the genomic bins (or probes) BED file.
- --fasta_file_path: Path to the genome FASTA file.
- --sample_labels_csv_path: Path to a CSV file containing sample labels, optional. Supplying this file allows you to differentiate line types for different labels in the `gc_bias_profile.png` output. The file should have two columns:
  - sample: Sample names matching the BAM file names.
  - \<label>: A column for your labels with "True" or "False" values.
- --out_dir: Path to the output directory.
- --at_anchor: GC percentile anchor for detecting AT bias. Should be > 0 and < 50. Default: 25
- --gc_anchor: GC percentile anchor for detecting GC bias. Should be > 50 and < 100. Default: 75
- --failure_fold_change: Relative coverage fold change failure threshold. Should be > 0. Default: 2
- --warning_fold_change: Relative coverage fold change warning threshold. Should be > 0 and less than failure_fold_change. Default: 1.5
- --failure_at: Coverage fold change failure threshold at the AT anchor. Should be > 0. Default: 1.5
- --failure_gc: Coverage fold change failure threshold at the GC anchor. Should be > 0. Default: 1.5
- --draw_trend: Boolean parameter to determine whether to generate trend visualization. Default: false
- --show_sample_names: Boolean parameter to determine whether to sample names in trend visualization. Default: true

## Memory Requirements
The bedtools_coverage process in panelGC is configured to use a maximum of 4 forks (panelGC.nf, line 130), as it typically requires around ~15GB per fork. It's important to note that users with less than 100GB of memory may need to decrease the number of forks, while those with more than 100GB can consider increasing it for potentially better performance.

To modify the fork settings, adjust the maxForks parameter in the bedtools_coverage process according to your system's memory capacity.

## Output
Files returned by panelGC:
1. gc_bias_loess_regression.tsv:
Records LOESS depth per GC percentile per sample.
2. gc_bias_loess_classification.tsv:
Records b<sub>25</sub>, b<sub>75</sub> and b<sub>75/25</sub> scores, and bias classification per sample.
3. gc_bias_profile.png:
The GC bias profile plot. Use the `--sample_labels_csv_path` argument to assign labels and differentiate line types for samples.
4. gc_bias_trend.png (optional):
The trend visualization of bias scores. Turned off by default.

## Demo Data

For demo data and test run instructions, please see [demo data documentation](demo_data/).

## Supplementary Scripts

### Simulate Paired End Sequencing Reads with GC/AT Bias

The `simulate_gc_bias.py` script simulates reads with GC/AT bias from a given reference genome and probe or genomic bins. 

For detailed usage instructions, see the [supplementary scripts documentation](supplementary_scripts/).

## Support
For support, questions, or contributions, please open an issue or a pull request in the repository or email easygsea@gmail.com

# Citation
Cheng X, Goktas MT, Williamson LM, Krzywinski M, Mulder DT, Swanson L, Slind J, Sihvonen J, Chow C, Carr A, Bosdet I, Tucker T, Young S, Moore R, Mungall KL, Yip S, Jones SJM. Enhancing clinical genomic accuracy with panelGC: a novel metric and tool for quantifying and monitoring GC biases in hybridization capture panel sequencing. Briefings in Bioinformatics. 2024 Sep;25(5):bbae442.
https://doi.org/10.1093/bib/bbae442
