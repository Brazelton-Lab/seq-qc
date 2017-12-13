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

import argparse
import bz2
import gzip
import io
from seq_qc import seq_io
import sys
from time import time

__author__ = "Christopher Thornton"
__license__ = 'GPLv2'
__maintainer__ = 'Christopher Thornton'
__status__ = "Production"
__version__ = "0.2.0"


def do_nothing(*args):
    pass


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fhandle', 
        metavar='in1.fast<q|a>', 
        action=seq_io.Open,
        mode='rb',
        help="input reads in fastq or fasta format. Can be a file containing "
             "either single-end or forward/interleaved reads if reads are "
             "paired-end [required]")
    input_arg = parser.add_mutually_exclusive_group(required=False)
    input_arg.add_argument('--interleaved',
        action='store_true',
        help="input is interleaved paired-end reads")
    input_arg.add_argument('-r', '--reverse',
        dest='rhandle',
        metavar='in2.fast<q|a>', 
        action=seq_io.Open,
        mode='rb',
        help="input reverse reads")
    parser.add_argument('-f', '--format',
        metavar='FORMAT',
        dest='format',
        default='fastq',
        choices=['fasta', 'fastq'],
        help="sequence file format. Can be fasta or fastq. [default: fastq]")
    parser.add_argument('-b', '--barcodes', 
        metavar='FILE',
        action=seq_io.Open,
        mode='r',
        help="file containing sample names mapped to the appropriate barcode "
             "sequences, in tab-separated format, with sample names in the "
             "first column. If this argument is unused, and the argument "
             "--force is given,the output files will be named for the "
             "barcode sequence found in the fasta\q file.")
    parser.add_argument('-s', '--suffix', 
        metavar='STR',
        type=str,
        help="string to append to the end of the file name. The default is to "
             "append the file format (fastq or fasta) and the strand for PE "
             "data (forward, reverse, interleaved).")
    parser.add_argument('--force',
        action='store_true',
        help="create new file for every barcode found in input")
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

    # Track program run-time
    start_time = time()


    # Assign variables based on arguments supplied by the user
    suffix = args.suffix if args.suffix else args.format
    if args.gzip:
        compression = '.gz'
        algo = GzipFile
    elif args.bzip2:
        compression = '.bz2'
        algo = BZ2File
    else:
        compression = ''
        algo = io.open


    # Prepare the iterator based on dataset type
    iterator = seq_io.read_iterator(args.fhandle, args.rhandle, \
                                   args.interleaved, args.format)


    # Store list of user-supplied barcodes
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


    # Demultiplex reads
    outfiles = {}
    unknowns = 0
    for processed_total, record in enumerate(iterator):
        # Prep output dependant on whether paired or unpaired
        try:
            tag = record.forward.description.split(':')[-1]
            ident = record.forward.id
            outf = record.forward.write()
            outr = record.reverse.write()
        except AttributeError:
            tag = record.description.split(':')[-1]
            ident = record.id
            outf = record.write()
            outr = None

        if (not tag.isalpha()) or (len(tag) != 6):
            seq_io.print_error("error: the format of the sequence headers is "
                               "incompatible with this method. Demultiplexing "
                               "these reads will require a different method "
                               "to be used instead")

        # Verify sequence tag in list of provided barcodes
        try:
            name = tags[tag]
        except KeyError:
            unknowns += 1
            if args.force:
                name = str(tag)
            else:
                print("warning: barcode {0} from sequence {1} doesn't correspond "
                      "to any of the barcodes provided. Use --force to write "
                      "sequences anyway".format(tag, ident), file=sys.stderr)
                continue

        # Write record to appropriate output file
        try:
            outfiles[name][0](outf)
            outfiles[name][1](outr)

        except KeyError:
            # Barcode not encountered previously, open new file for writes
            if args.rhandle:
                handle1 = io.TextIOWrapper(algo("{0}.forward.{1}{2}"\
                    .format(name, suffix, compression), mode='wb'))
                handle2 = io.TextIOWrapper(algo("{0}.reverse.{1}{2}"\
                    .format(name, suffix, compression), mode='wb'))
                write1, write2 = handle1.write, handle2.write
            elif args.interleaved:
                handle1 = io.TextIOWrapper(algo("{0}.interleaved.{1}{2}"\
                    .format(name, suffix, compression), mode='wb'))
                write1 = write2 = handle1.write
            else:
                handle1 = io.TextIOWrapper(algo("{0}.{1}{2}".format(name, \
                    suffix, compression), mode='wb'))
                write1 = handle1.write
                write2 = do_nothing

            outfiles[name] = (write1, write2)

            # Should be safe to write now
            outfiles[name][0](outf)
            outfiles[name][1](outr)


    # Verify input file non-empty
    try:
        processed_total += 1
    except UnboundLocalError:
        seq_io.print_error("error: no sequences were found to process")


    # Calculate and print output statistics
    num_parts = len(outfiles)
    print("\nRecords processed:\t{!s}\nBarcode partitions created:\t{!s}\n"
          "Sequences with unknown Barcode:\t{!s}\n".format(processed_total, \
          num_parts, unknowns), file=sys.stderr)


    # Calculate and print program run-time
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("It took {:.2e} minutes to process {!s} records\n"\
          .format(total_time, processed_total), file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
