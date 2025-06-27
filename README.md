# panelGC

## Introduction
panelGC effectively quantifies and monitors GC bias in hybridization capture sequencing, identifying and flagging potential sample and/or procedural anomalies.

## Supported Architectures
amd64, arm64v8, ppc64le, s390x ([more info](https://github.com/docker-library/official-images#architectures-other-than-amd64))

## Dependencies
- Apptainer (tested with 1.3.1)
- Nextflow (tested with 24.10.3)

Containerized dependencies (no installation required):
bedtools 2.30.0 HTSlib 1.22 SAMtools 1.22
r-base 4.3.2 argparser 0.7.1 tidyverse 2.0.0 data.table 1.14.8

## Installation
- Install Apptainer: \
Refer to official [Apptainer installation guide](https://apptainer.org/docs/user/main/quick_start.html)
- Install Nextflow: \
Refer to official [Nextflow installation guide](https://www.nextflow.io/docs/latest/getstarted.html)

## Usage
### Pull panelGC repository:
```bash
git clone https://github.com/easygsea/panelGC.git
```
### Run panelGC:
```bash
nextflow run /path/to/panelGC/panelGC.nf \
  --bam_directory_path /workspace/bam_files/ \
  --bed_file_path /workspace/bed_file.bed \
  --fasta_file_path /workspace/fasta_file.fa \
  --out_dir /workspace/output_directory/
```

**Note:** `nextflow` provides additional options to help optimize the deployment of panelGC for your specific use. For more details, please refer to the official guides: [nextflow CLI reference](https://www.nextflow.io/docs/latest/reference/cli.html#options).

### Parameters
- --bam_directory_path: Path to the directory containing alignment BAM files. Indices are preferred but not mandatory. Symlinks to the BAM and index files are valid.
- --bed_file_path: Path to the genomic bins (or probes) BED file.
- --fasta_file_path: Path to the genome FASTA file.
- --sample_labels_csv_path: Path to a CSV file containing sample labels, optional. Supplying this file allows you to differentiate line types for different labels in the `gc_bias_profile.png` output. The file should have two columns:
  - sample: Sample names matching the BAM file names.
  - \<label>: A column for your labels with "True" or "False" values.
- --out_dir: Path to the output directory.
- --window_size: Window size (bp) for calculating GC content. Regions in the BED file smaller than window_size will be skipped. Use 0 to calculate GC content for each entire region in the BED file. Default: 100
- --at_anchor: GC percentile anchor for detecting AT bias. Should be > 0 and < 50. Default: 25
- --gc_anchor: GC percentile anchor for detecting GC bias. Should be > 50 and < 100. Default: 75
- --failure_fold_change: Relative coverage fold change failure threshold. Should be > 0. Default: 2
- --warning_fold_change: Relative coverage fold change warning threshold. Should be > 0 and less than failure_fold_change. Default: 1.5
- --failure_at: Coverage fold change failure threshold at the AT anchor. Should be > 0. Default: 1.5
- --failure_gc: Coverage fold change failure threshold at the GC anchor. Should be > 0. Default: 1.5
- --y_lim: y-axis minimum and maximum for the GC bias profile plot. Comma-separated string of two numbers where the first number is less than the second e.g. "0,1". Default: "auto", which means the y-axis will be automatically determined by the data.
- --draw_trend: Boolean parameter to determine whether to generate trend visualization. Default: false
- --show_sample_names: Boolean parameter to determine whether to sample names in trend visualization. Default: true
- --draw_per_base_coverage: Boolean parameter to determine whether to draw per-base coverage plot. Default: true
- --publish_per_base_coverage: Boolean parameter to determine whether to publish the per-base coverage file for each sample in the output directory. Default: true
- --publish_gc_content_summary: Boolean parameter to determine whether to publish the GC content summary file in the output directory. Default: true
- --publish_bam_files: Boolean parameter to determine whether to publish converted BAM files in the output directory when input files are in CRAM format. Default: false

### Retry mechanism
The workflow uses a retry mechanism to handle transient failures. By default, processes will retry up to 10 times with an exponential backoff strategy. This helps ensure robust execution in environments with unstable network connections or resource constraints, particularly when processing CRAM files.

To modify retry behavior, you can adjust the following parameters in `nextflow.config`:
- `maxRetries`: Maximum number of retry attempts (default: 10)
- `errorStrategy`: Strategy for handling errors (default: 'retry')

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
