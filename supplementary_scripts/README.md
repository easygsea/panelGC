# Supplementary Scripts

## Simulate Paired End Sequencing Reads with GC/AT Bias

This script simulates paired end reads with GC/AT bias from a given reference genome and probe or genomic bins.

### Installation

#### Install Python 3.11.7

Download and install Python 3.11.7 from the official Python website: [https://www.python.org/downloads/release/python-3117/](https://www.python.org/downloads/release/python-3117/)

#### Set Up the Environment

To set up the environment, follow these steps:
1. Create a virtual environment
Follow the instructions from the official Python documentation: [Create virtual environments](https://docs.python.org/3/library/venv.html#creating-virtual-environments)

2. Activate the virtual environment:
Follow the instructions from the official Python documentation: [Activate virtual environments](https://docs.python.org/3/library/venv.html#how-venvs-work)

3. Install the dependencies:
```bash
pip install -r requirements.txt
```

### Usage

To run the script with parameters, follow the instructions from the official Python documentation to [activate your configured virtual environment](https://docs.python.org/3/library/venv.html#how-venvs-work), then run the following command:
```bash
python simulate_gc_bias.py path/to/reference.fasta path/to/probes.bed --read_length 150 --probe_length 120 --coverage 10 --bias_type GC --bias_degree 2.0
```

#### Parameters
- fasta_file: Input reference FASTA file
- bed_file: Input BED file with probe regions
- --read_length: Length of the reads to simulate (default: 150)
- --probe_length: Length of the probes (default: 120)
- --coverage: Coverage to simulate (default: 10)
- --bias_type: Type of bias ("GC" or "AT", default: "GC")
- --bias_degree: Degree of to apply (default: 2.0)
