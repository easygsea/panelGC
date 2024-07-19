# panelGC

## Introduction
panelGC effectively quantifies and monitors GC bias in hybridization capture sequencing, identifying and flagging potential sample and/or procedural anomalies.

## Supported Architectures
amd64, arm64v8, ppc64le, s390x ([more info](https://github.com/docker-library/official-images#architectures-other-than-amd64))

## Dependencies
- Nextflow (tested with 21.10.6)
- Singularity (tested with 2.5.1)

Containerized dependencies (no installation required):
bedtools 2.30.0
r-base 4.3.2 argparser 0.7.1 BiocManager 1.30.22 GenomicRanges 1.54.1 rtracklayer 1.62.0 tidyverse 2.0.0

## Installation
- Install Nextflow (if not already installed): \
Refer to official [Nextflow installation guide](https://www.nextflow.io/docs/latest/install.html)
- Install Singularity (if not included with Nextflow): \
Refer to official [Singularity installation guide](https://docs.sylabs.io/guides/latest/user-guide/quick_start.html)
- Clone the panelGC repository:
```bash
git clone https://github.com/easygsea/panelGC.git
```

## Usage
Run panelGC with the following command:
```bash
nextflow run panelGC.nf \
  --bam_directory_path /path/to/bam_files/ \
  --bed_file_path /path/to/bed_file.bed \
  --fasta_file_path /path/to/fasta_file.fa \
  --out_dir /path/to/output_directory/
```
Replace the paths with your actual data directories and file paths. Adjust [parameters](https://github.com/easygsea/panelGC/tree/main?tab=readme-ov-file#parameters) as needed.

Alternatively, configure parameters in a JSON or YML file (see official [Nextflow parameters guide](https://training.nextflow.io/basic_training/config/#parameters)) and excecute the following command:
```bash
nextflow run panelGC.nf -params-file params.json
```

**Note:** If you have limited space in your working directory or prefer to store the Singularity container in a different location, you can set the 'NXF_SINGULARITY_CACHEDIR' environment variable. This variable allows you to specify a custom path for storing Singularity images. To use this feature, export the 'NXF_SINGULARITY_CACHEDIR' variable with your desired path before running the panelGC command. For example:

```bash
export NXF_SINGULARITY_CACHEDIR=/path/to/singularity_cache/
```

### Parameters
- --bam_directory_path: Path to the directory containing alignment BAM files. Indices are preferred but not mandatory. Symlinks to the BAM and index files are valid.
- --bed_file_path: Path to the genomic bins (or probes) BED file.
- --fasta_file_path: Path to the genome FASTA file.
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
The bedtools_coverage process in panelGC is configured to use a maximum of 4 forks, as it typically requires around ~15GB per fork. It's important to note that users with less than 100GB of memory may need to decrease the number of forks, while those with more than 100GB can consider increasing it for potentially better performance.

To modify the fork settings, adjust the maxForks parameter in the bedtools_coverage process according to your system's memory capacity.

## Output
Files returned by panelGC:
1. gc_bias_loess_regression.tsv:
Records LOESS depth per GC percentile per sample.
2. gc_bias_loess_classification.tsv:
Records b<sub>25</sub>, b<sub>75</sub> and b<sub>75/25</sub> scores, and bias classification per sample.
3. gc_bias_profile.png:
The GC bias profile plot.
4. gc_bias_trend.png (optional):
The trend visualization of bias scores. Turned off by default.

## Demo Data

For demo data and test run nstructions, please see [demo data documentation](demo_data/).

## Supplementary Scripts

### Simulate Paired End Sequencing Reads with GC/AT Bias

The `simulate_gc_bias.py` script simulates reads with GC/AT bias from a given reference genome and probe or genomic bins. 

For detailed usage instructions, see the [supplementary scripts documentation](supplementary_scripts/).

## Support
For support, questions, or contributions, please open an issue or a pull request in the repository or email easygsea@gmail.com

# Citation
Cheng X, Goktas MT, Williamson LM, Krzywinski M, Mulder DT, Swanson L, Slind J, Sihvonen J, Chow C, Carr A, Bosdet I, Tucker T, Young S, Moore R, Mungall KL, Yip S, Jones SJM. Enhancing clinical genomic accuracy with panelGC: a novel metric and tool for quantifying and monitoring GC biases in hybridization capture panel sequencing. 2024.
