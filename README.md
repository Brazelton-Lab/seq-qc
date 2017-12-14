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

* arandomness
* bio_utils

## Installation

pip install seq_qc

## Usage

filter_replicates takes as input fastq or fasta files. Paired reads can either 
be split in separate files or can be in a single, interleaved file. 
filter_replicates currently supports searching for exact, 5'-prefix, and 
reverse-complement replicates.

### Examples

    filter_replicates --prefix --rev-comp --reverse input_reverse.fastq \
        input_forward.fastq

    filter_replicates --prefix --rev-comp -o output_forward.fastq.gz -v \
        output_reverse.fastq.gz --reverse input_reverse.fastq.gz \
        input_forward.fastq.gz

    filter_replicates --interleaved -o output_interleaved.fasta \
        --format fasta --log output.log input_interleaved.fasta

    filter_replicates -o output_singles.fastq input_singles.fastq

qtrim takes only fastq files as input. It can perform a variety of trimming 
steps - including trimming low quality bases from the start and end of a 
read, trimming the read after the position of the first ambiguous base, and
trimming a read using a sliding-window approach. It also supports cutting the 
read to a desired length by removing bases from the start or end of the read.

### Examples

    qtrim -1 -o \
    output_forward.fastq.gz --qual-type phred33 --sliding-window 10:20 -v \
        out_reverse.fastq.gz -s out_singles.fastq.gz --reverse \
        in_reverse.fastq.gz in_forward.fastq.gz

    qtrim --threads 4 --interleaved --leading 20 --trailing 20 --trunc-n \
        --min-len 60 -o out_interleaved.fastq in_interleaved.fastq

demultiplex_headers will split sequences into separate files by the barcode 
sequence found in the sequence headers. Headers must use the Casava 1.8 format
for demultiplexing to work. An optional barcodes file can be provided that maps
desired output filenames to the barcodes expected in the input Fastq/a file.
See demultiplex_headers --help for formatting requirements.

### Examples

    demultiplex_headers --barcodes barcode_map.tsv --distance 2 --hist \
        barcodes.hist --reverse in_reverse.fastq.gz input_forward.fastq.gz

    demultiplex_headers --bzip2 --interleaved input_interleaved.fastq.gz
