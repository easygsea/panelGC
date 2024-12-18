# Demo Data

This folder contains instructions and parameters file for a demo run to help you quickly try out panelGC.

## Data Preparation

Please follow these steps to download demo files into this demo_data/ folder:

1. Download hg19 genome fasta file from UCSC (https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/) and decompress as necessary. Example command:
```bash
wget --timestamping 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz'
gzip -d hg19.fa.gz
```
2. Download all 11 files in [panelGC demo alignment BAM files and probe BED file](https://dx.doi.org/10.6084/m9.figshare.26232365) from figshare and move them into this demo_data/ folder.
```bash
wget https://figshare.com/ndownloader/articles/26232365/versions/1
unzip 1
rm 1
```

After preparation, your demo_data/ folder should contain the following 13 files:
```bash
demo_data
├── demo_params.json
├── demo_probes.bed
├── hg19.fa
├── simulated_at_bias_high.bam
├── simulated_at_bias_high.bam.bai
├── simulated_at_bias_medium.bam
├── simulated_at_bias_medium.bam.bai
├── simulated_gc_bias_high.bam
├── simulated_gc_bias_high.bam.bai
├── simulated_gc_bias_medium.bam
├── simulated_gc_bias_medium.bam.bai
├── simulated_no_bias.bam
└── simulated_no_bias.bam.bai
```

## Run panelGC on demo data

Make sure you have properly [installed dependencies](https://github.com/easygsea/panelGC/tree/main?tab=readme-ov-file#installation). Go to your cloned panelGC parental folder which looks like this:

```bash
panelGC
├── bin
├── demo_data
├── nextflow.config
├── panelGC.nf
└── supplementary_scripts
```

Then execute:
```bash
singularity exec \
  -B $(pwd):/workspace panelgc_latest.sif \
  nextflow run panelGC.nf \
  -params-file demo_data/demo_params.json
```

## Example panelGC output from demo run

If executed successfully, you will notice a new demo_output/ folder:

```bash
demo_output
├── gc_bias_curves.png
├── gc_bias_loess_classification.tsv
├── gc_bias_loess_regression.tsv
└── gc_bias_trend.png

```

You may compare the output against [panelGC demo output](https://dx.doi.org/10.6084/m9.figshare.26232485) available on figureshare.