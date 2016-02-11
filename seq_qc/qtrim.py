#! /usr/bin/env python
"""
Trim and filter sequences using phred quality score information and length \
thresholds
"""

from __future__ import print_function
from __future__ import division

__author__ = "Christopher Thornton"
__date__ = "2016-01-18"
__version__ = "0.10.2"

import argparse
import seq_io
import sys
import trim

def apply_trimming(record, steps, qual_type, crop=None, headcrop=None, 
    trunc=False):
    seq = record['sequence']
    qual = record['quality']
    origlen = len(seq)

    if crop:
        seq, qual = (seq[:crop], qual[:crop])
    if headcrop:
        seq, qual = (seq[headcrop:], qual[headcrop:])

    for step, value in steps:
        start, end = step(qual, value, qual_type)
        seq, qual = (seq[start: end], qual[start: end])

    if trunc:
        start, end = trim.trunc_n(seq)
        seq, qual = (seq[start: end], qual[start: end])

    trimlen = len(seq)
    record['sequence'], record['quality'] = (seq, qual)
    return record, origlen, trimlen

def parse_sw_arg(argument):
    try:
        window, score = argument.split(':')
    except ValueError:
        seq_io.print_message("error: the input for -w/--window-size is "
            "formatted incorrectly. See --help for instructions", sys.stderr)
    else:
        if score.isdigit():
            score = int(score)
        else:
            seq_io.print_message("error: score threshold should be an integer "
                "value", sys.stderr)
        if window.isdigit():
            window = int(window)
        else:
            try:
                window = float(window)
            except ValueError:
                seq_io.print_message("error: window size should be either an "
                    "integer or a fraction", sys.stderr)

    return (window, score)

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-1', metavar='FILE', dest='f_file',
        required=True,
        help="input reads in fastq format. Can be a file containing either "
        "forward or interleaved reads if reads are paired-end [required]")
    input_arg = parser.add_mutually_exclusive_group(required=False)
    input_arg.add_argument('--interleaved',
        action='store_true',
        help="input is interleaved paired-end reads")
    input_arg.add_argument('-2', metavar='FILE', dest='r_file',
        help="input reverse reads in fastq format")
    parser.add_argument('-o', '--out', metavar='FILE', dest='out_f',
        required=True,
        type=seq_io.open_output,
        help="output trimmed reads [required]")
    parser.add_argument('-v', '--out-reverse', metavar='FILE', dest='out_r',
        type=seq_io.open_output,
        help="output trimmed reverse reads")
    parser.add_argument('-s', '--singles', metavar='FILE', dest='out_s',
        type=seq_io.open_output,
        help="output trimmed orphaned reads")
    parser.add_argument('-f', '--format', metavar='FORMAT', dest='out_format',
        default='fastq',
        choices=['fasta', 'fastq'],
        help="output files format (fastq or fasta) [default: fastq]")
    parser.add_argument('-l', '--log',
        type=seq_io.open_output,
        help="output log file to keep track of trimmed sequences")
    parser.add_argument('-q', '--qual-type', metavar='TYPE', dest='qual_type',
        default='phred33',
        choices=['phred33', 'phred64'],
        help="base quality encoding (phred33 or phred64) [default: phred33]")
    parser.add_argument('-m', '--min-len', metavar='LEN', dest='minlen',
        type=int,
        default=0,
        help="filter reads shorter than the threshold [default: 0]")
    trim_args = parser.add_argument_group('trimming options')
    trim_args.add_argument('-O', '--trim-order', metavar='ORDER',
        dest='trim_order',
        default='ltw',
        help="order of trimming steps [default: ltw (corresponds to leading, "
        "trailing, and sliding-window)]")
    trim_args.add_argument('-W', '--sliding-window', metavar='FRAME',
        dest='sw',
        type=parse_sw_arg,
        help="trim both 5' and 3' ends of a read using a sliding window "
        "approach. Input should be of the form 'window_size:qual_threshold', "
        "where 'qual_threshold' is an integer between 0 and 42 and "
        "'window_size' can either be length in bases or fraction of the total "
        "read length")
    trim_args.add_argument('-H', '--headcrop', metavar='BASES',
        type=int,
        help="remove exactly the number of bases specified from the start of "
        "the read")
    trim_args.add_argument('-C', '--crop', metavar='SIZE',
        type=int,
        help="trim to the specified size by removing bases from the end of "
        "the read")
    trim_args.add_argument('-L', '--leading', metavar='SCORE', 
        dest='lead_score',
        type=int,
        help="trim by removing low quality bases from the start of the read")
    trim_args.add_argument('-T', '--trailing', metavar='SCORE', 
        dest='trail_score',
        type=int,
        help="trim by removing low quality bases from the end of the read")
    trim_args.add_argument('--trunc-n', dest='trunc_n',
        action='store_true',
        help="truncate sequence at position of first ambiguous base")
    args = parser.parse_args()
    all_args = sys.argv[1:]

    iterator = seq_io.get_iterator(args.f_file, args.r_file, args.interleaved)

    if not args.out_s and (args.interleaved or args.r_file):
        parser.error("argument -s/--singles is required when arguments -2 or "
            "--interleaved are used")

    if args.r_file and not args.out_r:
        parser.error("argument -v/--out-reverse is required when argument -2 "
            "is used")

    seq_io.start_message('qtrim', all_args, __version__)

    trim_tasks = {'l': (trim.trim_leading, args.lead_score), 
        't': (trim.trim_trailing, args.trail_score), 
        'w': (trim.adaptive_trim, args.sw)}

    trim_steps = []
    for task in args.trim_order:
        value = trim_tasks[task][-1]
        if value:
            trim_steps.append(trim_tasks[task])
    if len(trim_steps) < 1 and not (args.crop or args.headcrop):
        seq_io.print_message("error: no trimming steps were applied", 
            sys.stderr)

    if args.out_format == 'fasta':
        writer = seq_io.fasta_writer
    else:
        writer = seq_io.fastq_writer

    paired = True if (args.interleaved or args.r_file) else False

    if paired:
        if args.interleaved:
            args.out_r = args.out_f

        seq_io.logger(args.log, "Record\tForward length\tForward trimmed "
            "length\tReverse length\tReverse trimmed length\n")

        num_pairs = pairs_passed = discarded_pairs = fsingles = rsingles = 0
        for forward, reverse in iterator:
            num_pairs += 1

            forward, flen, ftrim = apply_trimming(forward, trim_steps, 
                args.qual_type, args.crop, args.headcrop, args.trunc_n)
            reverse, rlen, rtrim = apply_trimming(reverse, trim_steps, 
                args.qual_type, args.crop, args.headcrop, args.trunc_n)

            # both good
            if ftrim >= args.minlen and rtrim >= args.minlen:
                pairs_passed += 1
                writer(args.out_f, forward)
                writer(args.out_r, reverse)
            # forward orphaned, reverse filtered
            elif ftrim >= args.minlen and rtrim < args.minlen:
                fsingles += 1
                writer(args.out_s, forward)
            # reverse orphaned, forward filtered
            elif ftrim < args.minlen and rtrim >= args.minlen:
                rsingles += 1
                writer(args.out_s, reverse)
            # both discarded
            else:
                discarded_pairs += 1

            seq_io.logger(args.log, "{}\t{}\t{}\t{}\t{}\n"
                .format(forward['identifier'], flen, ftrim, rlen, rtrim))

        total = num_pairs * 2
        passed = pairs_passed * 2 + fsingles + rsingles
        print("\nRecords processed:\t{!s} ({!s} pairs)\nPassed filtering:\t"
            "{!s} ({:.2%})\n  Paired reads kept:\t{!s} ({:.2%})\n  Forward "
            "only kept:\t{!s} ({:.2%})\n  Reverse only kept:\t{!s} ({:.2%})"
            "\nRead pairs discarded:\t{!s} ({:.2%})\n".format(total, num_pairs,
            passed, passed / total, pairs_passed, pairs_passed / num_pairs,
            fsingles, fsingles / total, rsingles, rsingles / total,
            discarded_pairs, discarded_pairs / num_pairs))

    else:
        seq_io.logger(args.log, "Record\tLength\tTrimmed length\n")

        total = discarded = 0
        for record in iterator:
            total += 1

            record, seqlen, trimlen = apply_trimming(record, trim_steps, 
                args.qual_type, args.crop, args.headcrop, args.trunc_n) 

            if trimlen >= args.minlen:
                writer(args.out_f, record)
            else:
                discarded += 1

            seq_io.logger(args.log, "{}\t{}\t{}\n".format(record['identifier'],
                seqlen, trimlen))

        passed = total - discarded
        print("\nRecords processed:\t{!s}\nPassed filtering:\t{!s} "
        "({:.2%})\nRecords discarded:\t{!s} ({:.2%})\n".format(total, passed,
        passed / total, discarded, discarded / total))

if __name__ == "__main__":
    main()
    sys.exit(0)
