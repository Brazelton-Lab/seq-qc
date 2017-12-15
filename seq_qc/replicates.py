#! /usr/bin/env python
"""
Remove replicate sequences from a FASTQ/A file. Supports searches for exact, 
5'-prefix, and reverse-complement exact and 5'-prefix replicates.
 
For split read pairs:
    filter_replicates [flags] -o out.forward [-v out.reverse] --r in.reverse \
        in.forward

For interleaved read pairs:
    filter_replicates [flags] [-o out.interleaved] in.interleaved

Supported file formats are FASTQ and FASTA. Compression using gzip and bzip2 
algorithms is automatically detected for input files. To compress the output 
files, add the appropriate file extension to the file names (.gz, .bz2). For 
interleaved or single-end reads, use /dev/stdin to indicate that input is from  
standard input (stdin). Similarly, leaving out '--output' will result in the 
output being sent to standard output (stdout).
"""

from __future__ import division
from __future__ import print_function

import argparse
from array import array
import hashlib
from seq_qc.seq_io import print_error, program_info, Open, read_iterator
from seq_qc.structure import reverse_complement, reverse_complement_paired
import sys
from time import time
import zlib

__author__ = "Christopher Thornton"
__license__ = 'GPLv2'
__maintainer__ = 'Christopher Thornton'
__status__ = "Production"
__version__ = "1.4.0"


def duplicate_status(seq1, seq2):
    """
    Determine if one sequence is a duplicate of another, and if so, whether \
    the duplicate is an exact match to or a prefix of the template sequence.

    A status of zero means not a duplicate, one means the second sequence is \
    an exact duplicate of the first, two means the second sequence is a \
    prefix of the first, and three means that the first sequence is a prefix \
    of the second.
    """
    seq1_len, seq2_len= (len(seq1), len(seq2))

    if seq1_len == seq2_len:
        if seq1 == seq2:
            return 1

    elif seq1_len > seq2_len:
        if seq1[:seq2_len] == seq2:
            return 2

    elif seq1_len < seq2_len:
        if seq1 == seq2[:seq1_len]:
            return 3

    return 0


def split_by_length(sequence, length):
    return sequence[:length], sequence[length:]


def do_nothing(*args):
    pass


def self(item):
    return item


def main():
    parser = argparse.ArgumentParser(description=__doc__,        
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fhandle',
        metavar='in1.fast<q|a>', 
        type=str,
        action=Open,
        mode='rb',
        default=sys.stdin,
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
        action=Open,
        mode='rb',
        help="input reverse reads")
    parser.add_argument('-o', '--out',
        metavar='FILE',
        dest='out_f',
        type=str,
        action=Open,
        mode='wt',
        default=sys.stdout,
        help="output trimmed reads [default: stdout]")
    parser.add_argument('-v', '--out-reverse', 
        metavar='FILE', 
        dest='out_r',
        type=str,
        action=Open,
        mode='wt',
        help="output reverse reads")
    parser.add_argument('-f', '--format', 
        metavar='FORMAT',
        dest='format',
        default='fastq',
        choices=['fasta', 'fastq'],
        help="sequence file format. Can be fasta or fastq. [default: fastq]")
    parser.add_argument('-l', '--log', 
        type=str,
        action=Open,
        mode='wt',
        help="output log of replicate types")
    dup_args = parser.add_argument_group('replicate types')
    dup_args.add_argument('--prefix',
        action='store_true',
        help="replicate can be a 5' prefix of another read")
    dup_args.add_argument('--rev-comp',
        dest='rev_comp',
        action='store_true',
        help="replicate can be the reverse-complement of another read")
    parser.add_argument('--reduce-memory', 
        dest='mem_use',
        action='store_true',
        help="reduce the mount of memory that the program uses. This could "
             "result in a drastic increase in run time.")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()
    all_args = sys.argv[1:]

    program_info('filter_replicates', all_args, __version__)

    # Track program run-time
    start_time = time()


    # Assign variables based on arguments supplied by the user
    out_f = args.out_f.write
    logger = args.log.write if args.log else do_nothing
    logger("#Replicate\tTemplate\tType\n")
    compress = zlib.compress if args.mem_use else self
    decompress = zlib.decompress if args.mem_use else self
    out_format = ">{0} {1}\n{2}\n" if args.format == "fasta" else \
                 "@{0} {1}\n{2}\n+\n{3}\n"
    paired = True if (args.interleaved or args.rhandle) else False


    # Prepare the iterator based on dataset type
    iterator = read_iterator(args.fhandle, args.rhandle, args.interleaved, \
                             args.format)


    # Assign variables based on dataset type (paired or single-end)
    if paired:
        print("Processing input as paired-end reads\n", file=sys.stderr)
        out_r = out_f if not args.out_r else args.out_r.write
        rc = reverse_complement_paired

    else:
        print("Processing input as single-end reads\n", file=sys.stderr)
        out_r = do_nothing
        rc = reverse_complement


    # Iterate over the reads, storing only the unique records
    uniques = {}
    for records_total, entry in enumerate(iterator):
        try:
            header = (entry.forward.id, entry.reverse.id)
            fdesc, rdesc = (entry.forward.description, entry.reverse.description)
            fseq, rseq = (entry.forward.sequence, entry.reverse.sequence)
        except AttributeError:
            header = (entry.id, '')
            fdesc, rdesc = (entry.description, '')
            fseq, rseq = (entry.sequence, '')
            try:
                qual = compress(entry.quality)
            except AttributeError:
                # Must be fasta format
                qual = None
        else:
            try:
                qual = compress(entry.forward.quality + entry.reverse.quality)
            except AttributeError:
                qual = None

        flen, rlen = len(fseq), len(rseq)
        record = [i for i in (header, flen, fseq + rseq, qual) if i != None]

        # Use hash of full or prefixed sequence as a key for quick comparisons
        fsub, rsub = ((20, 20) if args.prefix else (flen, rlen))
        key = hashlib.md5((fseq[:fsub] + rseq[:rsub]).encode()).digest()


        # Search if replicate
        search_list = []
        try:
            search_list = uniques[key]
        except KeyError:
            # No match to the database found. Need to check the reverse 
            # complement if requested
            if args.rev_comp:
                try:
                    fseq_rc, rseq_rc = rc(fseq, rseq)
                except TypeError:
                    fseq_rc, rseq_rc = rc(fseq), ''

                rckey = hashlib.md5((fseq_rc[:fsub] + rseq_rc[:rsub])\
                                    .encode()).digest()

                try:
                    search_list = uniques[rckey]
                except KeyError:
                    # Not a replicate. Add to the database of uniques
                    uniques[key] = [record]
                    continue
                else:
                    duplicate_key = rckey
                    fquery, rquery = fseq_rc, rseq_rc
                    query_id = header[0]

            # Not a replicate. Add to the database of uniques
            else:
                uniques[key] = [record]
                continue
        else:
            duplicate_key = key
            fquery, rquery = fseq, rseq
            query_id = header[0]


        # Search through list of records with common key to see if the sequence
        # matches one that has been observed before
        duplicate = None
        for index, search_record in enumerate(search_list):
            # Get search sequences by splitting combined sequence on forward 
            # length
            search_id = search_record[0][0]
            fsearch, rsearch = split_by_length(search_record[2], \
                                               search_record[1])

            # Check replicate status of forward sequence
            fstatus = duplicate_status(fquery, fsearch)
            if fstatus:
                # Check reverse only if forward a duplicate
                rstatus = duplicate_status(rquery, rsearch)
                if rstatus:
                    # Query is an exact match to a DB record
                    if (fstatus == 1 and rstatus == 1):
                        duplicate_type = "exact"
                        duplicate = query_id
                        template = search_id
                        break

                    # Query is a prefix of a DB record
                    elif (fstatus == 1 and rstatus == 3) or \
                         (fstatus == 3 and rstatus == 1) or \
                         (fstatus == 3 and rstatus == 3):
                        duplicate_type = "prefix"
                        duplicate = query_id
                        template = search_id
                        break

                    # A DB record is a prefix of the query
                    elif (fstatus == 1 and rstatus == 2) or \
                         (fstatus == 2 and rstatus == 1) or \
                         (fstatus == 2 and rstatus == 2):
                        duplicate_type = "prefix"
                        duplicate = search_id
                        template = query_id
                        # Replace old DB record with new
                        uniques[duplicate_key][index] = record
                        break
 
        if duplicate:
            # Add rc to duplicate type if search_list from rev_comp
            duplicate_type = "rev-comp {}".format(duplicate_type) \
                             if key != duplicate_key else duplicate_type

            logger("{}\t{}\t{}\n".format(duplicate, template, duplicate_type))

        else:
            # record is definitely not a duplicate, so add to the list of 
            # unique sequences with a common key
            uniques[duplicate_key].append(record)


    # Make sure input file non-empty
    try:
        records_total += 1  #number records processed
    except UnboundLocalError:
        print_error("error: no sequences were found to process.")


    # Write unique records
    uniques_total = 0  #remaining records after dereplication
    for unique_key in uniques:
        for record in uniques[unique_key]:
            uniques_total += 1
            fheader, rheader = record[0]

            fseq, rseq = split_by_length(record[2], record[1])
            try:
                fqual, rqual = split_by_length(decompress(record[3]), record[1])
            except IndexError:
                fqual = rqual = None

            out_f(out_format.format(*tuple([i for i in (fheader, fdesc, fseq, \
                                    fqual) if i != None])))
            out_r(out_format.format(*tuple([i for i in (rheader, rdesc, rseq, \
                                    rqual) if i != None])))


    # Calculate and print output statistics
    replicates_total = records_total - uniques_total
    print("Records processed:\t{!s}\nUnique reads found:\t{!s} ({:.2%})\nReplicate reads found:\t{!s} "
        "({:.2%})\n".format(records_total, uniques_total, uniques_total / records_total, replicates_total, replicates_total / records_total), file=sys.stderr)


    # Calculate and print program run-time
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("It took {:.2e} minutes to process {!s} records\n"\
          .format(total_time, records_total), file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
