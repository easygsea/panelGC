# panelGC

## Introduction
panelGC effectively quantifies and monitors GC bias in hybridization capture sequencing, identifying and flagging potential sample and/or procedural anomalies.

## Dependencies
- Nextflow (tested with 21.10.6)
- Singularity (tested with 2.5.1)

## Installation
- Install Nextflow (if not already installed):
```bash
curl -s https://get.nextflow.io | bash
```
- Install Singularity (if not included with Nextflow): \
Refer to the official [Singularity installation guide](https://docs.sylabs.io/guides/latest/user-guide/quick_start.html).
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
Replace the paths with your actual data directories and file paths.

## Memory Requirements
The bedtools_coverage process in panelGC is configured to use a maximum of 4 forks, as it typically requires around ~15GB per fork. It's important to note that users with less than 100GB of memory may need to decrease the number of forks, while those with more than 100GB can consider increasing it for potentially better performance.

To modify the fork settings, adjust the maxForks parameter in the bedtools_coverage process according to your system's memory capacity.

### Parameters
- --bam_directory_path: Path to the directory containing alignment BAM files.
- --bed_file_path: Path to the genomic bins (or probes) BED file.
- --fasta_file_path: Path to the genome FASTA file.
- --out_dir: Path to the output directory.

## Support
For support, questions, or contributions, please open an issue or a pull request in the repository.

## Citation
If you use panelGC in your research, please cite:
