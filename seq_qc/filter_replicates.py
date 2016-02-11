#! /usr/bin/env python
"""
De-replicate paired-end sequencing reads. Can check for exact, 5'-prefix, 
and reverse-complement duplicates.
"""
from __future__ import print_function
from __future__ import division

__author__ = "Christopher Thornton"
__date__ = "2016-01-17"
__version__ = "0.13.1"

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
    duplicate, two means template is a duplicate, and three means query is a \
    prefix duplicate.
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

def replicate_status(record, position, key, db, reverse={}):
    """
    Check if record is a prefix or exact duplicate of another read in the \
    dataset. The function can also be used to check the prefix \
    reverse-complement of a read or read pair if set.

    Returns the ID of the replicate, the ID of the template, the replicate \
    type, and the appropriate key in the database.
    """
    ident, fseq, rseq = record

    if key in db:
        query_set = db[key]
        for search_pos, search_obj in query_set:
            search_id, fsearch, rsearch = search_obj.split(',')
            fstatus = compare_seqs(fseq, fsearch)
            # check forward read first. If it is a duplicate then check reverse
            if fstatus:
                rstatus = compare_seqs(rseq, rsearch)
                if rstatus:
                    if (fstatus == 1 and rstatus == 1):
                        return (ident, search_id, 'exact', position)
                    elif (fstatus == 1 and rstatus == 3) or \
                        (fstatus == 3 and rstatus == 1) or \
                        (fstatus == 3 and rstatus == 3):
                        return (ident, search_id, 'prefix', position)
                    elif (fstatus == 1 and rstatus == 2) or \
                        (fstatus == 2 and rstatus == 1) or \
                        (fstatus == 2 and rstatus == 2):
                        return (search_id, ident, 'prefix', search_pos)

    return (None, None, None, None)

def main():
    parser = argparse.ArgumentParser(description=__doc__,        
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-1', metavar='FILE', dest='f_file',
        required=True,
        help="input reads in fastq or fasta format. Can be a file containing "
        "either forward or interleaved reads if reads are paired-end "
        "[required]")
    input_arg = parser.add_mutually_exclusive_group(required=False)
    input_arg.add_argument('--interleaved',
        action='store_true',
        help="input is interleaved paired-end reads")
    input_arg.add_argument('-2', metavar='FILE', dest='r_file',
        help="input reverse reads in fastq or fasta format")
    parser.add_argument('-o', '--out', dest='out_f', metavar='FASTQ',
        type=seq_io.open_output,
        help="output reads in fastq or fasta format")
    parser.add_argument('-v', '--out-reverse', metavar='FILE', dest='out_r',
        type=seq_io.open_output,
        help="output reverse reads in  fastq or fasta")
    parser.add_argument('-f', '--format', metavar='FORMAT',
        dest='out_format',
        default='fastq',
        choices=['fasta', 'fastq'],
        help="output files format (fastq or fasta) [default: fastq]")
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
    args = parser.parse_args()
    all_args = sys.argv[1:]

    iterator = seq_io.get_iterator(args.f_file, args.r_file, args.interleaved)

    seq_io.start_message('filter_replicates', all_args, __version__)
    seq_io.logger(args.log, "Replicate\tTemplate\tType\n")

    if args.out_format == 'fasta':
        writer = seq_io.fasta_writer
    else:
        writer = seq_io.fastq_writer

    paired = True if (args.interleaved or args.r_file) else False

    seq_db = {}
    duplicates = array('i')
    i = 0
    if paired:
        for forward, reverse in iterator:
            ident = forward['identifier']
            fseq, rseq = (forward['sequence'], reverse['sequence'])

            fsubsize, rsubsize = ((20, 20) if args.prefix else 
                (len(fseq), len(rseq)))
            key = hashlib.md5(fseq[:fsubsize] + rseq[:rsubsize]).digest()

            dup_id, dup_temp, dup_type, dup_pos = replicate_status(
                (ident, fseq, rseq), i, key, seq_db, reverse=reverse)
            if dup_id:
                duplicates.append(dup_pos)
                seq_io.logger(args.log, "{}\t{}\t{}\n".format(dup_id, dup_temp,
                    dup_type))
                i += 1
                continue

            # not a duplicate and reverse-complement checking is set
            if args.rev_comp:
                f_rc, r_rc = pairs.reverse_complement_paired(fseq, rseq)
                rckey = hashlib.md5(f_rc[:fsubsize] + r_rc[:rsubsize]).digest()
                dup_id, dup_temp, dup_type, dup_pos = replicate_status(
                    (ident, f_rc, r_rc), i, rckey, seq_db, reverse=reverse)
                if dup_id:
                    duplicates.append(dup_pos)
                    dup_type = 'rev-comp ' + dup_type
                    seq_io.logger(args.log, "{}\t{}\t{}\n".format(dup_id, dup_temp,
                        dup_type))
                    i += 1
                    continue

            # record is definitely not a duplicate, so add to database
            try:
                seq_db[key].append((i, "{},{},{}".format(ident, fseq, rseq)))
            except KeyError:
                seq_db[key] = [(i, "{},{},{}".format(ident, fseq, rseq))]
            i += 1

    else:
        for record in iterator:
            ident = record['identifier']
            seq = record['sequence']

            subsize = 20 if args.prefix else len(seq)
            key = hashlib.md5(seq[:subsize]).digest()

            dup_id, dup_temp, dup_type, dup_pos = replicate_status(
                (ident, seq, ''), i, key, seq_db)
            if dup_id:
                duplicates.append(dup_pos)
                seq_io.logger(args.log, "{}\t{}\t{}\n".format(dup_id, dup_temp,
                    dup_type))
                i += 1
                continue

            # not a duplicate and reverse-complement checking is set
            if args.rev_comp:
                rc = reverse_complement(seq)
                rckey =  hashlib.md5(rc[:subsize]).digest()
                dup_id, dup_temp, dup_type, dup_pos = replicate_status(
                    (ident, rc, ''), i, rckey, seq_db)
                if dup_id:
                    duplicates.append(dup_pos)
                    dup_type = 'rev-comp ' + dup_type
                    seq_io.logger(args.log, "{}\t{}\t{}\n".format(dup_id, dup_temp,
                        dup_type))
                    i += 1
                    continue

            try:
                seq_db[key].append((i, "{},{},{}".format(ident, seq, '')))
            except KeyError:
                seq_db[key] = [(i, "{},{},{}".format(ident, seq, ''))]
            i += 1

    if args.out_f:
        iterator = seq_io.get_iterator(args.f_file, args.r_file,
            args.interleaved)
        duplicates = sorted(duplicates)
        index = 0
        try:
            dup_pos = duplicates[index]
        except IndexError:
            print("No duplicates found")
        else:
            if paired:
                if args.interleaved:
                    args.out_r = args.out_f

                for position, (forward, reverse) in enumerate(iterator):
                    if position != dup_pos:
                        writer(args.out_f, forward)
                        writer(args.out_r, reverse)
                    else:
                        index += 1
                        try:
                            dup_pos = duplicates[index]
                        except IndexError:
                            continue
            else:
                for position, record in enumerate(iterator):
                    if position != dup_pos:
                        writer(args.out_f, record)
                    else:
                        index += 1
                        try:
                            dup_pos = duplicates[index]
                        except IndexError:
                            continue

    num_dups = len(duplicates)
    print("\nRecords processed:\t{!s}\nReplicates found:\t{!s} "
        "({:.2%})\n".format(i, num_dups, num_dups / i))

if __name__ == "__main__":
    main()
    sys.exit(0)
