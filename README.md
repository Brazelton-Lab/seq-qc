# seq_qc - Python package with modules for preprocessing sequencing reads

## About

seq_qc is a python package for performing various quality control tasks on 
sequencing reads. Currently, seq_qc has three programs for this - a 
dereplicator of paired-end reads, a tool for performing quality-based
trimming, and a demultiplexer that uses the barcodes in sequence headers to
partition sequences.

## Requirements

Python 3.4+

Python Libraries:

* screed
* arandomness
* bio_utils

## Installation

pip install seq_qc

## Usage

filter_replicates takes as input fastq or fasta files. Paired reads can either 
be split in separate files or can be in a single, interleaved file. The types 
of replicates that filter_replicates can search for are exact, 5'-prefix, and 
reverse-complement replicates.

### Examples

    filter_replicates -1 input_forward.fastq -2 input_reverse.fastq --prefix \
    --rev-comp

    filter_replicates -1 input_forward.fastq.gz -2 input_reverse.fastq.gz -o \
    output_forward.fastq.gz -v output_reverse.fastq.gz --prefix --rev-comp

    filter_replicates -1 input_interleaved.fastq -o output_interleaved.fasta \
    --interleaved --format fasta --log output.log

    filter_replicates -1 input_singles.fastq -o output_singles.fastq

qtrim takes only fastq files as input. It can perform a variety of trimming 
steps - including trimming low quality bases from the start and end of a 
read, trimming the read after the position of the first ambiguous base, and
trimming a read using a sliding-window approach. It also supports cutting the 
read to a desired length by removing bases from the start or end of the read.

### Examples

    qtrim -1 input_forward.fastq.gz -2 input_reverse.fastq.gz -o \
    output_forward.fastq.gz -v output_reverse.fastq.gz -s \
    output_singles.fastq.gz --qual-type phred33 --sliding-window 10:20

    qtrim -1 input_interleaved.fastq -o output_interleaved.fastq --interleaved \
    --leading 20 --trailing 20 --trunc-n --min-len 60

### Examples

    demultiplex_headers --barcodes barcode_map.tsv input_forward.fastq.gz \
    input_reverse.fastq.gz

    demultiplex_headers --bzip2 --interleaved input_interleaved.fastq.gz
