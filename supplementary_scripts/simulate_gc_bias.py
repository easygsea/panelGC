#! /usr/bin/python

import argparse
import gzip
import random
from Bio import SeqIO
import numpy as np

def read_fasta(fasta_file):
    with open(fasta_file, 'r') as file:
        return SeqIO.to_dict(SeqIO.parse(file, 'fasta'))

def read_bed(bed_file):
    probes = []
    with open(bed_file, 'r') as file:
        for line in file:
            chr_name, start, end, gene_name = line.strip().split()
            start, end = int(start), int(end)
            probes.append((chr_name, start, end, gene_name))
    return probes

def gc_content(seq):
    return float(seq.count('G') + seq.count('C')) / len(seq)

def at_content(seq):
    return float(seq.count('A') + seq.count('T')) / len(seq)

def simulate_reads_with_bias(reference, probes, read_length, coverage, bias_type, bias_degree, probe_length):
    output_prefix = 'simulated_reads_' + bias_type.lower() + str(bias_degree)
    read_overhang = int((read_length - probe_length) / 2)

    with gzip.open(f'{output_prefix}_1.fq.gz', 'wt') as f1, gzip.open(f'{output_prefix}_2.fq.gz', 'wt') as f2:
        for chr_name, start, end, gene_name in probes:
            seq = reference[chr_name].seq[start - read_overhang:start + read_length + read_overhang].upper()

            seq_probe = seq[read_overhang:read_overhang + probe_length]
            if bias_type == 'GC':
                content = gc_content(seq_probe)
            elif bias_type == 'AT':
                content = at_content(seq_probe)
            else:
                raise ValueError('Invalid bias type. Use "GC" or "AT".')

            if content > 0.5:
                num_reads = coverage * (1 + random.gauss(1, 0.1) * (content - 0.5) * bias_degree)
            else:
                num_reads = coverage * (1 - random.gauss(1, 0.1) * (0.5 - content) * bias_degree)

            for i in range(int(num_reads)):
                read_start = random.randint(0, 29)
                read1 = seq[read_start:read_start + read_length]
                read2 = read1.reverse_complement()
                f1.write(f'@{gene_name}_{i}\n{read1}\n+\n{"~" * read_length}\n')
                f2.write(f'@{gene_name}_{i}\n{read2}\n+\n{"~" * read_length}\n')

def main():
    parser = argparse.ArgumentParser(description='Simulate reads with GC/AT bias.')
    parser.add_argument('fasta_file', type=str, help='Input reference FASTA file')
    parser.add_argument('bed_file', type=str, help='Input BED file with probe regions')
    parser.add_argument('--read_length', type=int, default=150, help='Length of the reads to simulate (default: 150)')
    parser.add_argument('--probe_length', type=int, default=120, help='Length of the probes (default: 120)')
    parser.add_argument('--coverage', type=int, default=10, help='Coverage to simulate (default: 10)')
    parser.add_argument('--bias_type', type=str, choices=['GC', 'AT'], default='GC', help='Type of bias ("GC" or "AT", default: "GC")')
    parser.add_argument('--bias_degree', type=float, default=2.0, help='Degree of bias to apply (default: 2.0)')

    args = parser.parse_args()

    reference = read_fasta(args.fasta_file)
    probes = read_bed(args.bed_file)
    simulate_reads_with_bias(reference, probes, args.read_length, args.coverage, args.bias_type, args.bias_degree, args.probe_length)

if __name__ == '__main__':
    main()