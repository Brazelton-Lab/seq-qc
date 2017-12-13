#! /usr/bin/env python
"""
Split sequences into separate files based on the barcode sequence found in
the sequence header. The sequence headers should be of the form 
'@seqid <strand>:N:0:<barcode>' (Casava 1.8).

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
from bz2 import BZ2File
from gzip import GzipFile
import io
from seq_qc import seq_io
import sys
from time import time

__author__ = "Christopher Thornton"
__license__ = 'GPLv2'
__maintainer__ = 'Christopher Thornton'
__status__ = "Production"
__version__ = "0.2.0"


class BarcodeEntry(object):
    """A simple class to store template barcodes

    Attributes:
            id (str): Barcode identifier

            sequence (str): Barcode sequence

            count (int): Number of times template barcode has been observed
    """
    def __init__(self, initval=0):
        """Initialize attributes to store Barcode entry data"""
        self.id = None
        self.sequence = None
        self.count = initval

    def increment(self):
        self.count += 1


def do_nothing(*args):
    pass


def hamming_distance(s1, s2):
    #Return Hamming distance between equal-length sequences
    if len(s1) != len(s2):
        raise ValueError("Template and sequence index barcodes must be equal "
                         "in length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def return_last(last):
    """Return last item in a list or tuple"""
    return last[-1]


def sort_by_last(tuple_list):
    """Sort list of tuples or lists by the value of their last item"""
    return sorted(tuple_list, reverse=False, key=return_last)


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
    parser.add_argument('-d', '--distance',
        type=int,
        default=0,
        help="hamming distance allowed between sequence barcodes in order to "
             "be placed into the same partition. Requires a barcodes file "
             "providing template barcode sequences.")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()
    all_args = sys.argv[1:]

    seq_io.program_info('demultiplex_headers', all_args, __version__)

    if args.distance and not args.barcodes:
        parser.error("error: argument -b/--barcodes must be used with "
                     "-d/--distance")

    # Track program run-time
    start_time = time()


    # Assign variables based on arguments supplied by the user
    if args.barcodes and args.distance > 0:
        outstats = "Records processed:\t{0}\nBarcode partitions created:\t{1}"\
                   "\nSequence barcodes with -\n  exact match to a template:"\
                   "\t{2}\n  one or more mismatchs:\t{3}\nSequences with "\
                   "unknown barcode:\t{4}\n"
        exact_total = mismatch_total = unknowns = 0
    elif args.barcodes and args.distance == 0:
        outstats = "Records processed:\t{0}\nBarcode partitions created:\t{1}"\
                   "\nSequences with unknown barcode:\t{2}\n"
        exact_total = mismatch_total = None
        unknowns = 0
    else:
        outstats = "Records processed:\t{0}\nBarcode partitions created:\t{1}"\
                   "\n"
        exact_total = mismatch_total = unknowns = None

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
        for line in args.barcodes:
            tag = BarcodeEntry()

            # Verify barcodes file correctly formatted
            try:
                name, sequence = line.strip().split('\t')
            except ValueError:
                seq_io.print_error("error: barcode mapping file does not "
                                   "appear to be formatted correctly")

            # Verify unique sample names
            if name in [i.id for i in tags.values()]:
                seq_io.print_error("error: the same sample name is used for "
                                   "more than one barcode sequence")

            tag.id = name
            tag.sequence = sequence

            tags[sequence] = tag


    # Demultiplex reads
    outfiles = {}
    for processed_total, record in enumerate(iterator):
        # Prepare output dependant on whether paired or unpaired
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


        # Find the template barcode with the smallest hamming distance to the 
        # record sequence barcode
        if args.distance:
            distances = sort_by_last([(i.sequence, hamming_distance(tag, \
                                     i.sequence)) for i in tags.values()])

            min_tag, min_dist = distances[0]

            if min_dist == 0:
                exact_total += 1
            else:
                mismatch_total += 1

            # Determine if more than one closest match
            if [i[1] for i in distances].count(min_dist) > 1:
                seq_io.print_warning("warning: barcode {0} in sequence {1} is "
                                     "equally similar to more than one "
                                     "template barcode. Unable to determine "
                                     "which partition to assign it to"\
                                     .format(tag, ident))
                continue
            else:
                if min_dist <= args.distance:
                    tag = min_tag
 

        # Verify sequence tag in list of provided barcodes
        if args.barcodes:
            try:
                file_prefix = tags[tag].id
            except KeyError:
                unknowns += 1
                if args.force:
                    file_prefix = str(tag)
                else:
                    seq_io.print_warning("warning: sequence barcode {0} does "
                                         "not correspond to any of the "
                                         "template barcodes provided. Use "
                                         "--force to write these records "
                                         "anyway".format(tag))
                    continue

        else:
            file_prefix = str(tag)


        # Write record to appropriate output file
        try:
            outfiles[file_prefix][0](outf)
            outfiles[file_prefix][1](outr)

        except KeyError:
            # Barcode not encountered previously, open new file for writes
            if args.rhandle:
                handle1 = io.TextIOWrapper(algo("{0}.forward.{1}{2}"\
                    .format(file_prefix, suffix, compression), mode='wb'))
                handle2 = io.TextIOWrapper(algo("{0}.reverse.{1}{2}"\
                    .format(file_prefix, suffix, compression), mode='wb'))
                write1, write2 = handle1.write, handle2.write
            elif args.interleaved:
                handle1 = io.TextIOWrapper(algo("{0}.interleaved.{1}{2}"\
                    .format(file_prefix, suffix, compression), mode='wb'))
                write1 = write2 = handle1.write
            else:
                handle1 = io.TextIOWrapper(algo("{0}.{1}{2}".format(file_prefix, \
                    suffix, compression), mode='wb'))
                write1 = handle1.write
                write2 = do_nothing

            outfiles[file_prefix] = (write1, write2)

            # Should be safe to write now
            outfiles[file_prefix][0](outf)
            outfiles[file_prefix][1](outr)


    # Verify input file non-empty
    try:
        processed_total += 1
    except UnboundLocalError:
        seq_io.print_error("error: no sequences were found to process")


    # Calculate and print output statistics
    partitions_total = len(outfiles)
    stats = [processed_total, partitions_total] + [i for i in \
             (exact_total, mismatch_total, unknowns) if i != None]
    print(outstats.format(*tuple(stats)), file=sys.stderr)


    # Calculate and print program run-time
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("It took {:.2e} minutes to process {!s} records\n"\
          .format(total_time, processed_total), file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
