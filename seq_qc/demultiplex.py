#! /usr/bin/env python
"""
Split sequences into separate files based on the barcode sequences found in
the headers. Headers should be formatted as Casava 1.8 (i.e 
'@seqid <strand>:N:0:<barcode>').

For single-end and interleaved reads:
    demultiplex_headers [options] input
 
For split paired-end reads:
    demultiplex_headers [options] in.forward in.reverse

Supported file formats are FASTQ and FASTA. Compression using gzip and bzip2 
algorithms is automatically detected for the input files through their file 
extension. For single-end or interleaved reads, use /dev/stdin to indicate 
that input is from standard input (stdin). Similarly, leaving out '--output' 
will result in the output being sent to standard output (stdout).
"""

from __future__ import print_function

import argparse
from bz2 import BZ2File
from gzip import GzipFile
import io
from seq_qc import seq_io
from statistics import mean, median, stdev
import sys
from time import time

__author__ = "Christopher Thornton"
__license__ = 'GPLv2'
__maintainer__ = 'Christopher Thornton'
__status__ = "Production"
__version__ = "1.0.1"


class BarcodeEntry(object):
    """A simple class to store template barcodes

    Attributes:
            id (str): Barcode identifier

            barcode (str): Barcode sequence and multiplex group, if applicable

            count (int): Number of times template barcode has been observed
    """
    def __init__(self, initval=0):
        """Initialize attributes to store Barcode entry data"""
        self.id = None
        self.barcode = None
        self.group = None
        self.count = initval

    def increment(self):
        self.count += 1

    def sequence(self):
        return self.barcode.split(':')[0]

    def run_info(self):
        try:
            return ':'.join(self.barcode.split(':')[1:])
        except IndexError:
            return ''


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
        help="sequence file format [default: fastq]. Available options are "
             "'fasta' or 'fastq'")
    parser.add_argument('-b', '--barcodes', 
        metavar='FILE',
        action=seq_io.Open,
        mode='r',
        help="file containing sample names mapped to template barcode "
             "sequences, in tab-separated format. The first column should "
             "contain sample names and the second column should contain the "
             "appropriate barcodes. An optional third column can be used to "
             "assign barcodes to 'multiplex groups'. The name of a multiplex "
             "group should have first the run id, then the flowcell lane, "
             "separated by a colon (e.g. 432:4). If this argument is unused "
             "and --force is provided instead, the output files will be named "
             "for the barcode and run information found in the headers")
    parser.add_argument('-s', '--suffix', 
        metavar='STR',
        type=str,
        help="string to append to the end of the file name. The default is to "
             "append the file format (fastq or fasta) and the strand for PE "
             "data (forward, reverse, interleaved)")
    parser.add_argument('-c', '--hist',
        metavar='FILE',
        action=seq_io.Open,
        mode='w',
        help="output histogram of barcode counts. This can be used for "
             "graphing the error distribution of a barcode sequence, for "
             "instance")
    parser.add_argument('--no-out',
        dest='no_out',
        action='store_true',
        help="do not write sequences to a file. Only output a histogram of "
             "barcode counts")
    parser.add_argument('--force',
        action='store_true',
        help="create new file for every barcode found in input")
    compress_arg = parser.add_mutually_exclusive_group(required=False)
    compress_arg.add_argument('--gzip',
        action='store_true',
        help="output files should be compressed using the gzip algorithm. The "
             "suffix '.gz'. will be appended to the file names")
    compress_arg.add_argument('--bzip2',
        action='store_true',
        help="output files should be compressed using the bzip2 algorithm. The "
             "suffix '.bz2' will be appended to the file names")
    parser.add_argument('-d', '--distance',
        type=int,
        default=0,
        help="maximum hamming distance allowed between sequence barcodes to "
             "be sorted into the same partition. Requires a barcodes file "
             "providing template barcode sequences")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()
    all_args = sys.argv[1:]

    seq_io.program_info('demultiplex_by_header', all_args, __version__)

    if args.distance and not args.barcodes:
        parser.error("error: argument -d/--distance cannot be used without "
                     "-b/--barcodes")

    if args.no_out and not args.hist:
        parser.error("error: argument --no-out cannot be used without "
                     "-c/--hist")

    # Track program run-time
    start_time = time()


    # Assign variables based on arguments supplied by the user
    if args.barcodes and args.distance > 0:
        outstats = "Records processed:\t{0}\nBarcode partitions created:\t{1}"\
                   "\nSequence barcodes with -\n  exact match to a template:"\
                   "\t{2}\n  one or more mismatchs:\t{3}\nSequences with "\
                   "unknown barcode:\t{4}\n"
    elif args.barcodes and args.distance == 0:
        outstats = "Records processed:\t\t{0}\nBarcode partitions created:\t{1}"\
                   "\nSequences with unknown barcode:\t{2}\n"
    else:
        outstats = "Records processed:\t\t{0}\nBarcode partitions created:\t{1}"\
                   "\n"

    suffix = args.suffix if args.suffix else args.format
    out_hist = args.hist.write if args.hist else do_nothing
    exact_total = mismatch_total = unknowns = 0

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


    # Store user-supplied barcodes
    template_barcodes = {}
    names = []
    is_grouped = False
    if args.barcodes:
        barcode_dists = []
        for barcodes_total, line in enumerate(args.barcodes):
            tag = BarcodeEntry()
            group = ''

            # Verify that barcodes file is correctly formatted
            try:
                # Has filename and barcode only
                name, barcode = line.strip().split('\t')
            except ValueError:
                try:
                    # Has filename, barcode, and multiplex group
                    name, barcode, group = line.strip().split('\t')
                    is_grouped = True
                except ValueError:
                    seq_io.print_error("error: barcode mapping file does not "
                                       "appear to be formatted correctly. See "
                                       "the help message for formatting "
                                       "requirements")

            tag.id = name
            tag.barcode = barcode
            tag.group = group

            # Add to dictionary of template barcodes
            try:
                dict_group = template_barcodes[group]
            except KeyError:
                template_barcodes[group] = {barcode: tag}
            else:
                # Verify unique sample names and calculate barcode distances
                for i in dict_group:
                    if name == dict_group[i].id:
                        seq_io.print_error("error: sample name '{}' is being "
                                           "used for more than one template "
                                           "barcode".format(name))

                    barcode_dists.append((group, hamming_distance(dict_group[i]\
                                         .barcode, barcode)))

                if barcode not in dict_group:
                    dict_group[barcode] = tag
                else:
                    seq_io.error("error: template barcode '{}' has been seen "
                                 "more than once in a single multiplex group"\
                                 .format(barcode))


        # Verify non-empty barcodes file
        try:
            barcodes_total += 1
        except UnboundLocalError:
            seq_io.print_error("error: no barcodes were found to process")


        # Template barcode statistics
        out_bstats = "Template barcodes found:\t{!s}\n  # multiplex groups:\t{!s}\n\n"

        num_tem_groups = len(template_barcodes)
        bstats = [barcodes_total, num_tem_groups]

        if is_grouped:
            out_bstats += "  Multiplex Group\t# Barcodes\tMin distance\tMax distance\n"
            for mul_group in template_barcodes:
                out_bstats += "  {!s}\t{!s}\t{!s}\t{!s}\n"

                group_dists = [return_last(i) for i in barcode_dists if \
                               i[0] == group]

                bstats += [mul_group, len(template_barcodes[mul_group]), min(group_dists), \
                           max(group_dists)]

        else:
            out_bstats += "Minimum Hamming distance between all barcodes:\t"\
                          "{!s}\nMaximum Hamming distance between all barcodes:"\
                          "\t{!s}\n"

            all_dists = [return_last(i) for i in barcode_dists]

            bstats += [min(all_dists), max(all_dists)]

        print(out_bstats.format(*tuple(bstats)), file=sys.stderr)


    # Demultiplex reads
    outfiles = {}
    sequence_barcodes = []
    for processed_total, record in enumerate(iterator):
        # Prepare output dependant on whether paired or unpaired
        try:
            seq_tag = record.forward.description.split(':')[-1]
            header = record.forward.id
            outf = record.forward.write()
            outr = record.reverse.write()
        except AttributeError:
            seq_tag = record.description.split(':')[-1]
            header = record.id
            outf = record.write()
            outr = None


        # Verify headers are formatted as Casava 1.8
        try:
            run_info = tuple(header.split(':')[1: 4])  #run id, flowcell, lane
        except IndexError:
            seq_io.print_error("error: the format of the sequence headers is "
                               "incompatible with this method. Demultiplexing "
                               "these reads will require a different method "
                               "to be used instead")

        if (not seq_tag.isalpha()) or (len(seq_tag) != 6):
            seq_io.print_error("error: the format of the sequence headers is "
                               "incompatible with this method. Demultiplexing "
                               "these reads will require a different method "
                               "to be used instead")


        file_prefix = "{0}_{1}_{2}_{3}".format(seq_tag, *run_info)

        # Increment barcode occurences
        index = "{0}:{1}:{2}:{3}".format(seq_tag, *run_info)
        already_seen = False
        for i in sequence_barcodes:
            if index == i.barcode:
                i.increment()
                already_seen = True

        if not already_seen:
            seq_entry = BarcodeEntry()
            seq_entry.barcode = index
            seq_entry.count = 1
            sequence_barcodes.append(seq_entry)


        # Map sequence barcodes to barcodes in the provided barcodes file
        if args.barcodes:
            # Only consider other barcodes within the same multiplex group
            if is_grouped:
                try:
                    relevant_barcodes = template_barcodes[":".join(run_info)]
                except KeyError:
                    seq_io.warning("warning: run information in sequence "
                                   "header '{}' doesn't match any multiplex "
                                   "group in the provided barcodes file\n"\
                                   .format(header))
                    continue

            else:
                # Consider all barcodes within file
                relevant_barcodes = template_barcodes['']


            # Find the template barcode with the smallest hamming distance to 
            # the record sequence barcode
            distances = sort_by_last([(i, hamming_distance(seq_tag, i)) \
                                     for i in relevant_barcodes])

            min_tag, min_dist = distances[0]  #minimum is the first element

            if min_dist == 0:
                exact_total += 1
            else:
                mismatch_total += 1

            # Determine if more than one closest match
            if [i[1] for i in distances].count(min_dist) > 1:
                seq_io.print_warning("warning: barcode {0} in sequence "
                                     "{1} is equally similar to more than "
                                     "one template barcode. Unable to "
                                     "determine which partition to assign "
                                     "it to".format(seq_tag, header))
                continue
            else:
                # Assign to template partition if within threshold distance
                if min_dist <= args.distance and args.distance:
                    seq_tag = min_tag
 

            # Verify sequence tag is in the list of provided barcodes
            try:
                file_prefix = relevant_barcodes[seq_tag].id
            except KeyError:
                unknowns += 1
                if not args.force:
                    seq_io.print_warning("warning: barcode {0} in sequence {1}"
                                         " does not correspond to any of the "
                                         "template barcodes provided. The "
                                         "template ({2}) with the fewest "
                                         "mismatches is {3} nucleotides "
                                         "different. Use --force to write "
                                         "these records anyway".format(seq_tag,\
                                         header, min_tag, min_dist))
                    continue

        
        # Write record to appropriate output file
        if not args.no_out:
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
                    handle1 = io.TextIOWrapper(algo("{0}.{1}{2}".format(
                        file_prefix, suffix, compression), mode='wb'))
                    write1 = handle1.write
                    write2 = do_nothing

                outfiles[file_prefix] = (write1, write2)

                # Should be safe to write now
                outfiles[file_prefix][0](outf)
                outfiles[file_prefix][1](outr)


    # Write output histogram and sequence barcode statistics
    houtstats = "Sequence barcodes found:\t{!s}\n  Mean abundance:\t"\
                "{:.2}\n  Median abundance:\t{:.2}\n  Abundance SD:\t"\
                "\t{:.2}\n\n  Barcode\tRun Information\tAbundance\n"

    num_seq_tags = len(sequence_barcodes)
    houtstats += "  {!s}\t{!s}\t{!s}\n" * num_seq_tags

    barcode_abundances = []
    hstats_extra = []
    for seq_bar in sequence_barcodes:
        abund = seq_bar.count
        barcode_abundances.append(abund)
        hstats_extra += [seq_bar.sequence(),seq_bar.run_info(), abund]


    b_mean, b_median, b_sd = (mean(barcode_abundances), \
                              median(barcode_abundances), \
                              stdev(barcode_abundances))

    hstats = [num_seq_tags, b_mean, b_median, b_sd] + hstats_extra
    print(houtstats.format(*tuple(hstats)), file=sys.stderr)

    if args.hist:
        out_hist("#Total:  {}\n#Mean:   {:.2}\n#Median: {:.2}\n#STDev:  "
                 "{:.2}\n".format(num_seq_tags, b_mean, b_median, b_sd))

        for abundance in sorted(set(barcode_abundances)):
            counts = barcode_abundances.count(abundance)
            out_hist("{0}\t{1}\n".format(abundance, counts))


    # Verify input file non-empty
    try:
        processed_total += 1
    except UnboundLocalError:
        seq_io.print_error("error: no sequences were found to process")


    # Calculate and print output statistics
    if not args.no_out:
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
