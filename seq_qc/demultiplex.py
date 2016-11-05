#! /usr/bin/env python
"""
Split sequences into separate files based on the barcode sequence found in
the sequence header. The sequence headers should be of the form 
'@seqid <strand>:N:0:<barcode>'.

For single-end and interleaved reads:
    demultiplex_headers [options] input
 
For split paired-end reads:
    demultiplex_headers [options] in.forward in.reverse

Supported file formats are FASTQ and FASTA. Compression using gzip and bzip2 
algorithms is automatically detected for the input files. For single-end or
interleaved reads, use '-' to indicate that input should be taken from standard
input (stdin).
"""

from __future__ import print_function

__author__ = "Christopher Thornton"
__date__ = "2016-10-29"
__version__ = "0.1.0"

import argparse
import seq_io
import sys

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('f_file', metavar='in1.fast<q|a>', 
        help="input reads in fastq or fasta format. Can be a file containing "
        "either single-end or forward/interleaved reads if reads are "
        "paired-end [required]")
    input_arg = parser.add_mutually_exclusive_group(required=False)
    input_arg.add_argument('--interleaved',
        action='store_true',
        help="input is interleaved paired-end reads")
    input_arg.add_argument('r_file', metavar='in2.fast<q|a>', nargs='?',
        help="input reverse reads")
    parser.add_argument('-b', '--barcodes', metavar='FILE',
        type=seq_io.open_input,
        help="file containing sample names mapped to the appropriate barcode"
        "sequences, in tab-separated format, with sample names in the first "
        "column. If this argument is unused, the output files will be named "
        "for each barcode sequence found in the fasta\q file.")
    parser.add_argument('-s', '--suffix', metavar='STR',
        type=str,
        help="string to append to the end of the file name. The default is to "
        "append the file format (fastq or fasta) and the strand for PE data "
        "(forward, reverse, interleaved).")
    parser.add_argument('-f', '--out-format', metavar='FORMAT',
        dest='out_format',
        default='fastq',
        choices=['fasta', 'fastq'],
        help="output file format. Can be fasta or fastq. [default: fastq]")
    compress_arg = parser.add_mutually_exclusive_group(required=False)
    compress_arg.add_argument('--gzip',
        action='store_true',
        help="output files should be compressed using the gzip algorithm. The "
        "suffix '.gz'. will be appended to the file names.")
    compress_arg.add_argument('--bzip2',
        action='store_true',
        help="output files should be compressed using the bzip2 algorithm. The "
        "suffix '.bz2' will be appended to the file names.")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)

    args = parser.parse_args()
    all_args = sys.argv[1:]

    seq_io.program_info('demultiplex_headers', all_args, __version__)

    f_file = sys.stdin if args.f_file == '-' else args.f_file
    iterator = seq_io.get_iterator(f_file, args.r_file, args.interleaved)

    writer = seq_io.fasta_writer if args.out_format == 'fasta' else \
            seq_io.fastq_writer

    if args.gzip:
        compression = '.gz'
    elif args.bzip2:
        compression = '.bz2'
    else:
        compression = ''

    suffix = args.suffix if args.suffix else args.out_format

    tags = {}
    if args.barcodes:
        names = []
        for line in args.barcodes:
            try:
                name, tag = line.strip().split('\t')
            except ValueError:
                seq_io.print_error("error: barcode mapping file does not "
                    "appear to be formatted correctly")
            if name in names:
                seq_io.print_error("error: the same sample name is used for "
                    "more than one barcode sequence")
            else:
                names.append(name)

            tags[tag] = name

    outfiles = {}
    for i, (forward, reverse) in enumerate(iterator):
        tag = forward['description'].split(':')[-1]
        if (not tag.isalpha()) or (len(tag) != 6):
            seq_io.print_error("error: unable to determine the format of the "
                "sequence headers")

        try:
            name = tags[tag]
        except KeyError:
            name = str(tag)

        try:
            writer(outfiles[name][0], forward)
            writer(outfiles[name][1], reverse)
        except KeyError:
            if args.r_file:
                handle1 = seq_io.open_output("{}.forward.{}{}".format(name,
                    suffix, compression))
                handle2 = seq_io.open_output("{}.reverse.{}{}".format(name,
                    suffix, compression))
            elif args.interleaved:
                handle1 = seq_io.open_output("{}.interleaved.{}{}".format(name,
                     suffix, compression))
                handle2 = handle1
            else:
                handle1 = seq_io.open_output("{}.{}{}".format(name, suffix,
                    compression))
                handle2 = ''
            outfiles[name] = [handle1, handle2]

            writer(handle1, forward)
            writer(handle2, reverse)

    i += 1
    num_parts = len(outfiles)
    print("\nRecords processed:\t{!s}\nNumber of partitions:\t{!s}\n".format(i, 
        num_parts), file=sys.stderr)

if __name__ == "__main__":
    main()
    sys.exit(0)
