#! /usr/bin/env python
"""
De-replicate paired-end sequencing reads. Can check for exact, 
5'-prefix, and reverse-complement replicates.
 
For split paired-end reads:
    filter_replicates [flags] -o out.forward -v out.reverse in.forward in.reverse

For interleaved paired-end reads:
    filter_replicates [flags] [-o out.interleaved] in.interleaved

Supported file formats are FASTQ and FASTA. Compression using gzip and bzip2 
algorithms is automatically detected for input files. To compress output, add
the appropriate file extension to the file names (.gz, .bz2). For interleaved
reads, use the file name '-' to indicate that input should be taken from 
standard input (stdin). Similarly, leaving out the -o argument will cause 
output to be sent to standard output (stdout).
"""
from __future__ import print_function
from __future__ import division

__author__ = "Christopher Thornton"
__date__ = "2016-04-25"
__version__ = "1.0.15"

import argparse
from array import array
import hashlib
import pairs
import seq_io
import sys
from screed.dna import reverse_complement

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

def main():
    parser = argparse.ArgumentParser(description=__doc__,        
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('f_file', metavar='in1.fast<q|a>', 
        help="input forward or interleaved reads [required]")
    input_arg = parser.add_mutually_exclusive_group(required=True)
    input_arg.add_argument('--interleaved',
        action='store_true',
        help="input is interleaved paired-end reads")
    input_arg.add_argument('r_file', metavar='in2.fast<q|a>', nargs='?',
        help="input reverse reads")
    parser.add_argument('-o', '--out', dest='out_f', metavar='FILE',
        type=seq_io.open_output, default=sys.stdout,
        help="output reads")
    parser.add_argument('-v', '--out-reverse', metavar='FILE', dest='out_r',
        type=seq_io.open_output,
        help="output reverse reads")
    parser.add_argument('-f', '--out-format', metavar='FORMAT',
        dest='out_format',
        default='fastq',
        choices=['fasta', 'fastq'],
        help="output file format. Can be fasta or fastq. [default: fastq]")
    parser.add_argument('-l', '--log', metavar='LOG',
        type=seq_io.open_output,
        help="output log file to keep track of replicates")
    dup_args = parser.add_argument_group('replicate types')
    dup_args.add_argument('--prefix',
        action='store_true',
        help="replicate can be a 5' prefix of another read")
    dup_args.add_argument('--rev-comp', dest='rev_comp',
        action='store_true',
        help="replicate can be the reverse-complement of another read")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()
    all_args = sys.argv[1:]

    if args.r_file and not args.out_r:
        parser.error("argument -v/--out-reverse is required when an input "
            "reverse file is provided")

    f_file = sys.stdin if args.f_file == '-' else args.f_file
    iterator = seq_io.get_iterator(f_file, args.r_file, args.interleaved)

    seq_io.logger(args.log, "Replicate\tTemplate\tType\n")

    if args.out_format == 'fasta':
        writer = seq_io.fasta_writer
    else:
        writer = seq_io.fastq_writer

    seq_db = {}
    uniques = {}
    for i, (forward, reverse) in enumerate(iterator):
        ident = forward['identifier']
        fdesc, rdesc = (forward['description'], reverse['description'])
        fseq, rseq = (forward['sequence'], reverse['sequence'])
        fqual, rqual = (forward['quality'], reverse['quality'])

        flen, rlen = len(fseq), len(rseq)

        uniques[i] = (fseq + rseq, flen, fqual + rqual, ident)

        fsubsize, rsubsize = ((20, 20) if args.prefix else (flen, rlen))
        key = hashlib.md5(fseq[:fsubsize] + rseq[:rsubsize]).digest()

        dup_pos, temp_pos, dup_type = replicate_status(i, key, uniques, seq_db)

        # match to database found, so delete id from database of uniques
        if dup_pos:
            seq_io.logger(args.log, "{}\t{}\t{}\n".format(uniques[dup_pos][3], 
                uniques[temp_pos][3], dup_type))
            try:
                del uniques[dup_pos]
            except KeyError:
                print("input file has more than one sequence with the same "
                    "identifier", sys.stderr)
                sys.exit(1)
            continue

        # sequence is unique, so check reverse-complement if set
        if args.rev_comp:
            f_rc, r_rc = pairs.reverse_complement_paired(fseq, rseq)
            rckey = hashlib.md5(f_rc[:fsubsize] + r_rc[:rsubsize]).digest()
            dup_pos, temp_pos, dup_type = replicate_status(i, rckey,  uniques,
                seq_db)
            if dup_pos:
                dup_type = 'rev-comp ' + dup_type
                seq_io.logger(args.log, "{}\t{}\t{}\n".format(
                    uniques[dup_pos][3], uniques[temp_pos][3], dup_type))
                try:
                    del uniques[dup_pos]
                except KeyError:
                    print("input file has more than one sequence with the same "
                        "identifier", sys.stderr)
                    sys.exit(1)
                continue

        # record is definitely not a duplicate, so add to database of ids to 
        # check a match for
        try:
            seq_db[key].append(i)
        except KeyError:
            seq_db[key] = [i]

    try:
        i += 1
    except UnboundLocalError:
        seq_io.print_error("error: no sequences were found to process")

    if args.interleaved:
        args.out_r = args.out_f

    for j, i in enumerate(sorted(uniques.keys())):
        record = uniques[i]
        ident = record[3]
        fseq, rseq = split_by_length(record[0], record[1])
        fqual, rqual = split_by_length(record[2], record[1])
        writer(args.out_f, {'identifier': ident, 'description': fdesc, 'sequence': fseq, 'quality': fqual})
        writer(args.out_r, {'identifier': ident, 'description': rdesc, 'sequence': rseq, 'quality': rqual})

    j += 1

    seq_io.program_info('filter_replicates', all_args, __version__)
    num_reps = i - j
    print("\nRead Pairs processed:\t{!s}\nReplicates found:\t{!s} "
        "({:.2%})\n".format(i, num_reps, num_reps / i), file=sys.stderr)

if __name__ == "__main__":
    main()
    sys.exit(0)
