#! /usr/bin/env python
"""
Interleave paired-end sequencing reads
"""

from __future__ import print_function

__author__ = "Christopher Thornton"
__date__ = "2016-01-16"
__version__ = "0.2.3"

import argparse
from itertools import izip
import seq_io
import sys

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-1', metavar='FILE', dest='f_iter',
        required=True,
        type=seq_io.open_input,
        help="input forward reads [required]")
    parser.add_argument('-2', metavar='FILE', dest='r_iter',
        required=True,
        type=seq_io.open_input,
        help="input reverse reads [required]")
    parser.add_argument('-o', '--out', metavar='FILE', dest='out_f',
        required=True,
        type=seq_io.open_output,
        help="output interleaved reads [required]")
    parser.add_argument('-f', '--format', metavar='FORMAT', dest='out_format',
        default='fastq',
        choices=['fasta', 'fastq'],
        help="output file format (fastq or fasta) [default: fastq]")
    args = parser.parse_args()
    all_args = sys.argv[1:]

    if args.out_format == 'fasta':
        writer = seq_io.fasta_writer
    else:
        writer = seq_io.fastq_writer

    seq_io.start_message('interleave_pairs', all_args, __version__)

    num_pairs = 0
    iterator = izip(args.f_iter, args.r_iter)
    for forward, reverse in iterator:
        num_pairs += 1
        writer(args.out_f, forward)
        writer(args.out_f, reverse)

    print("\nRecords processed:\t{!s} ({!s} pairs)\n".format(num_pairs * 2, 
        num_pairs))

if __name__ == "__main__":
    main()
    sys.exit(0)
