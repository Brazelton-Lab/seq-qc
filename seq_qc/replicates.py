#! /usr/bin/env python
"""
De-replicate paired-end sequencing reads. Can check if a read is an exact, 
5'-prefix, or reverse-complement duplicate of another read in the dataset.
 
For split read pairs:
    filter_replicates [flags] -o out.forward -v out.reverse in.forward in.reverse

For interleaved read pairs:
    filter_replicates [flags] [-o out.interleaved] in.interleaved

Supported file formats are FASTQ and FASTA. Compression using gzip and bzip2 
algorithms is automatically detected for input files. To compress output, add
the appropriate file extension to the file names (.gz, .bz2). For interleaved
reads, use the file name '-' to indicate that input should be taken from 
standard input (stdin). Similarly, leaving out the -o argument will cause 
output to be sent to standard output (stdout).
"""

from __future__ import division
from __future__ import print_function

import argparse
from array import array
import hashlib
from seq_qc import pairs, seq_io
from screed.dna import reverse_complement
import sys
from time import time
import zlib

__author__ = "Christopher Thornton"
__license__ = 'GPLv2'
__maintainer__ = 'Christopher Thornton'
__status__ = "Production"
__version__ = "1.4.0"


def compare_seqs(query, template):
    """
    Return the replicate status of a search.

    A status of zero means not a duplicate, one means query is exactly \
    duplicate, two means template is a prefix duplicate, and three means \
    query is a prefix duplicate.
    """
    query_len, temp_len= (len(query), len(template))

    if query_len == temp_len:
        if query == template:
            return 1
    elif query_len > temp_len:
        if query[:temp_len] == template:
            return 2
    elif query_len < temp_len:
        if query == template[:query_len]:
            return 3

    return 0


def split_by_length(sequence, length):
    return sequence[:length], sequence[length:]


def replicate_status(query_position, key, unique_db, search_db):
    """
    Check if record is a prefix or exact duplicate of another read in the \
    dataset. The function can also be used to check the prefix \
    reverse-complement of a read or read pair if set.

    Returns the ID of the replicate, the ID of the template, and what type \
    of replicate was found.
    """
    query_record = unique_db[query_position]
    fquery, rquery = split_by_length(query_record[0], query_record[1])

    if key in search_db:
        for search_position in search_db[key]:
            try:
                search_record = unique_db[search_position]
            # id deleted from uniques already, so skip
            except KeyError:
                continue
            fsearch, rsearch = split_by_length(search_record[0], search_record[1])
            fstatus = compare_seqs(fquery, fsearch)
            # check forward read first. If it is a duplicate then check reverse
            if fstatus:
                rstatus = compare_seqs(rquery, rsearch)
                if rstatus:
                    if (fstatus == 1 and rstatus == 1):
                        return (query_position, search_position, 'exact')
                    elif (fstatus == 1 and rstatus == 3) or \
                        (fstatus == 3 and rstatus == 1) or \
                        (fstatus == 3 and rstatus == 3):
                        return (query_position, search_position, 'prefix')
                    elif (fstatus == 1 and rstatus == 2) or \
                        (fstatus == 2 and rstatus == 1) or \
                        (fstatus == 2 and rstatus == 2):
                        return (search_position, query_position, 'prefix')

    return (None, None, None)


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
        action=seq_io.Open,
        mode='rb',
        default=sys.stdin,
        help="input reads in fastq or fasta format. Can be a file containing "
             "either single-end or forward/interleaved reads if reads are "
             "paired-end [required]")
    input_arg = parser.add_mutually_exclusive_group(required=True)
    input_arg.add_argument('--interleaved',
        action='store_true',
        help="input is interleaved paired-end reads")
    input_arg.add_argument('-r', '--reverse',
        dest='rhandle', 
        metavar='in2.fast<q|a>', 
        action=seq_io.Open,
        mode='rb',
        help="input reverse reads")
    parser.add_argument('-o', '--out',
        metavar='FILE',
        dest='out_f',
        type=str,
        action=seq_io.Open,
        mode='wt',
        default=sys.stdout,
        help="output trimmed reads [default: stdout]")
    parser.add_argument('-v', '--out-reverse', 
        metavar='FILE', 
        dest='out_r',
        type=str,
        action=seq_io.Open,
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
        action=seq_io.Open,
        mode='wt',
        help="output log file to keep track of replicates")
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

    seq_io.program_info('filter_replicates', all_args, __version__)

    # Track program run-time
    start_time = time()


    # Assign variables based on arguments supplied by the user
    out_f = args.out_f.write
    logger = args.log.write if args.log else do_nothing
    compress = zlib.compress if args.mem_use else self
    decompress = zlib.decompress if args.mem_use else self
    out_r = out_f if not args.out_r else args.out_r.write
    out_format = ">{0} {1}\n{2}\n" if args.format == "fasta" else \
                 "@{0} {1}\n{2}\n+\n{3}\n"

    # Prepare the iterator based on dataset type
    iterator = seq_io.read_iterator(args.fhandle, args.rhandle, \
                                    args.interleaved, args.format)


    logger("Replicate\tTemplate\tType\n")


    # Iterate over records, storing unique records in uniques
    seq_db = {}
    uniques = {}
    for i, record in enumerate(iterator):
        ident = (record.forward.id, record.reverse.id)
        fdesc, rdesc = record.forward.description, record.reverse.description
        fseq, rseq = (record.forward.sequence, record.reverse.sequence)
        fqual, rqual = (record.forward.quality, record.reverse.quality)

        flen, rlen = len(fseq), len(rseq)

        uniques[i] = (fseq + rseq, flen, compress(fqual + rqual), ident)

        fsubsize, rsubsize = ((20, 20) if args.prefix else (flen, rlen))
        key = hashlib.md5((fseq[:fsubsize] + rseq[:rsubsize]).encode()).digest()

        dup_pos, temp_pos, dup_type = replicate_status(i, key, uniques, seq_db)

        # match to database found, so delete id from database of uniques
        if dup_pos:
            logger("{}\t{}\t{}\n".format(uniques[dup_pos][3], \
                   uniques[temp_pos][3], dup_type))
            try:
                del uniques[dup_pos]
            except KeyError:
                seq_io.print_error("error: input file has more than one "
                    "sequence with the same identifier")
                sys.exit(1)
            continue

        # sequence is unique, so check reverse-complement if set
        if args.rev_comp:
            f_rc, r_rc = pairs.reverse_complement_paired(fseq, rseq)
            rckey = hashlib.md5((f_rc[:fsubsize] + r_rc[:rsubsize]).encode())\
                .digest()
            dup_pos, temp_pos, dup_type = replicate_status(i, rckey,  uniques,
                seq_db)
            if dup_pos:
                dup_type = 'rev-comp ' + dup_type
                logger("{}\t{}\t{}\n".format(uniques[dup_pos][3], \
                       uniques[temp_pos][3], dup_type))
                try:
                    del uniques[dup_pos]
                except KeyError:
                    seq_io.print_error("error: input file has more than one "
                        "sequence with the same identifier")
                continue

        # record is definitely not a duplicate, so add to database of ids to 
        # check a match for
        try:
            seq_db[key].append(i)
        except KeyError:
            seq_db[key] = [i]


    # Make sure input file non-empty
    try:
        i += 1  #number records processed
    except UnboundLocalError:
        seq_io.print_error("error: no sequences were found to process.")


    # Write unique records
    for j, index in enumerate(sorted(uniques.keys())):
        record = uniques[index]
        fident, rident = record[3]

        fseq, rseq = split_by_length(record[0], record[1])
        fqual, rqual = split_by_length(decompress(record[2]), record[1])


        out_f(out_format.format(fident, fdesc, fseq, fqual))
        out_r(out_format.format(rident, rdesc, rseq, rqual))

    j += 1  #number remaining records after dereplication


    # Calculate and print output statistics
    num_reps = i - j
    print("\nRead Pairs processed:\t{!s}\nReplicates found:\t{!s} "
        "({:.2%})\n".format(i, num_reps, num_reps / i), file=sys.stderr)


    # Calculate and print program run-time
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("It took {:.2e} minutes to process {!s} records\n"\
          .format(total_time, i), file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
