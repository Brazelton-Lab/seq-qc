#! /usr/bin/env python
"""
Standard sequence quality control tool that can be used for base cropping, 
trimming sequences by quality score, and filtering reads that fail to meet
a minimum length threshold post cropping and trimming.

For single-end and interleaved reads:
    qtrim [options] [-o output] input
 
For split paired-end reads:
    qtrim [option] -o out.forward -v out.reverse -s out.singles in.forward 
        in.reverse

Input files must be in FASTQ format. Output can be in either FASTQ or FASTA 
format. Compression using gzip and bzip2 algorithms is automatically detected 
for input files. To compress output, add the appropriate file extension to the 
file names (.gz, .bz2). For single-end or interleaved reads, use '-' to 
indicate that input should be taken from standard input (stdin). Similarly, 
leaving out the -o argument will cause output to be sent to standard output 
(stdout).
"""

from __future__ import print_function
from __future__ import division

__author__ = "Christopher Thornton"
__date__ = "2016-11-06"
__version__ = "1.5.1"

import argparse
import seq_io
import sys
import trim

def apply_trimming(record, steps, qual_type, start=0, end=0, trunc=False):
    qual = record['quality']
    seq = record['sequence']
    slen = len(seq)

    scores = trim.translate_quality(qual, qual_type)
    for step, value in steps:
        newstart, newend = step(scores[start: slen - end], value)
        start += newstart
        end += newend

    last = slen - end
    if trunc:
        try:
            nstart = seq[start: last].index('N')
        except ValueError:
            nstart = last
        seq, qual = seq[start: nstart], qual[start: nstart]
    else:
        seq, qual = seq[start: last], qual[start: last]

    record['sequence'], record['quality'] = seq, qual
    return record

def parse_sw_arg(argument):
    try:
        window, score = argument.split(':')
    except ValueError:
        seq_io.print_error("error: the input for -w/--window-size is "
            "formatted incorrectly. See --help for instructions")
    else:
        if score.isdigit():
            score = int(score)
        else:
            seq_io.print_error("error: score threshold should be an integer "
                "value")
        if window.isdigit():
            window = int(window)
        else:
            try:
                window = float(window)
            except ValueError:
                seq_io.print_error("error: window size should be either an "
                    "integer or a fraction")

    return (window, score)

def get_list(argument):
    try:
        argument = [abs(int(i.lstrip())) for i in argument.split(",")]
    except ValueError:
        seq_io.print_error("error: input to -c/--crop, -d/--headcrop, and "
            "-m/--min-len must be in the form INT or INT,INT")
    arglen = len(argument)
    if arglen < 1 or arglen > 2:
        seq_io.print_error("error: one or two integer values should be "
            "provided with -c/--crop, -h/--headcrop, and -m/--min-len")
    return argument

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('f_file', 
        metavar='in1.fastq',
        help="input reads in fastq format. Can be a file containing either "
        "single-end or forward/interleaved reads if reads are paired-end "
        "[required]")
    input_arg = parser.add_mutually_exclusive_group(required=False)
    input_arg.add_argument('--interleaved',
        action='store_true',
        help="input is interleaved paired-end reads")
    input_arg.add_argument('--force',
        action='store_true',
        help="force process as single-end reads even if input is interleaved "
        "paired-end reads")
    input_arg.add_argument('r_file', 
        metavar='in2.fastq', 
        nargs='?',
        help="input reverse reads in fastq format")
    parser.add_argument('-o', '--out', 
        metavar='FILE', 
        dest='out_f',
        type=seq_io.open_output, 
        default=sys.stdout,
        help="output trimmed reads [required]")
    output_arg = parser.add_mutually_exclusive_group(required=False)
    output_arg.add_argument('-v', '--out-reverse', 
        metavar='FILE', 
        dest='out_r',
        type=seq_io.open_output,
        help="output trimmed reverse reads")
    output_arg.add_argument('--out-interleaved', 
        dest='out_interleaved',
        action='store_true',
        help="output interleaved paired-end reads, even if input is split")
    parser.add_argument('-s', '--singles', 
        metavar='FILE', 
        dest='out_s',
        type=seq_io.open_output,
        help="output trimmed orphaned reads")
    parser.add_argument('-f', '--out-format', 
        metavar='FORMAT', 
        dest='out_format', 
        choices=['fasta', 'fastq'],
        default='fastq',
        help="output files format (fastq or fasta) [default: fastq]")
    parser.add_argument('-l', '--log',
        type=seq_io.open_output,
        help="output log file to keep track of trimmed sequences")
    parser.add_argument('-q', '--qual-type', 
        metavar='TYPE', 
        dest='qual_type',
        type=int, 
        choices=[33, 64],
        default=33,
        help="ASCII base quality score encoding [default: 33]. Options are "
            "33 (for phred33) or 64 (for phred64)")
    parser.add_argument('-m', '--min-len', 
        metavar='LEN [,LEN]', 
        dest='minlen',
        type=get_list, 
        default='0',
        help="filter reads shorter than the minimum length threshold [default:"
             " off]. Different values can be provided for the forward and "
             "reverse reads, respectively, by separating them with a comma "
             "(e.g. 80,60), or a single value can be provided for both")
    trim_args = parser.add_argument_group('trimming options')
    trim_args.add_argument('-O', '--trim-order', 
        metavar='ORDER',
        dest='trim_order',
        default='ltw',
        help="order of trimming steps [default: ltw (corresponds to leading, "
        "trailing, and sliding-window)]")
    trim_args.add_argument('-W', '--sliding-window', 
        metavar='FRAME',
        dest='sw',
        type=parse_sw_arg,
        help="trim both 5' and 3' ends of a read using a sliding window "
        "approach. Input should be of the form 'window_size:qual_threshold', "
        "where 'qual_threshold' is an integer between 0 and 42 and "
        "'window_size' can either be length in bases or fraction of the total "
        "read length")
    trim_args.add_argument('-H', '--headcrop', 
        metavar='INT [,INT]',
        type=get_list, 
        default='0',
        help="remove exactly the number of bases specified from the start of "
        "the reads. Different values can be provided for the forward and "
        "reverse reads, respectively, by separating them with a comma "
        "(e.g. 2,0), or a single value can be provided for both")
    trim_args.add_argument('-C', '--crop', 
        metavar='INT [,INT]',
        type=get_list, 
        default='0',
        help="remove exactly the number of bases specified from the end of "
        "the reads. Different values can be provided for the forward and "
        "reverse reads, respectively, by separating them with a comma "
        "(e.g. 2,0), or a single value can be provided for both")
    trim_args.add_argument('-L', '--leading', 
        metavar='SCORE', 
        dest='lead_score',
        type=int,
        help="trim by removing low quality bases from the start of the read")
    trim_args.add_argument('-T', '--trailing', 
        metavar='SCORE', 
        dest='trail_score',
        type=int,
        help="trim by removing low quality bases from the end of the read")
    trim_args.add_argument('--trunc-n', 
        dest='trunc_n',
        action='store_true',
        help="truncate sequence at the position of the first ambiguous base")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()
    all_args = sys.argv[1:]

    seq_io.program_info('qtrim', all_args, __version__)

    try:
        fcrop, rcrop = args.crop
    except ValueError:
        fcrop = rcrop = args.crop[0]
    try:
        fheadcrop, rheadcrop = args.headcrop
    except ValueError:
        fheadcrop = rheadcrop = args.headcrop[0]
    try:
        fminlen, rminlen = args.minlen
    except ValueError:
        fminlen = rminlen = args.minlen[0]

    f_file = sys.stdin if args.f_file == '-' else args.f_file
    out_f = args.out_f
    iterator = seq_io.get_iterator(f_file, args.r_file, args.interleaved)

    if args.r_file and not (args.out_r or args.out_interleaved):
        parser.error("one of -v/--out-reverse or --out-interleaved is required "
            "when the argument -r/--reverse is used")

    trim_tasks = {'l': (trim.trim_leading, args.lead_score), 
        't': (trim.trim_trailing, args.trail_score), 
        'w': (trim.adaptive_trim, args.sw)}

    trim_steps = []
    for task in args.trim_order:
        value = trim_tasks[task][-1]
        if value:
            trim_steps.append(trim_tasks[task])
    if len(trim_steps) < 1 and not (args.crop or args.headcrop):
        seq_io.print_error("error: no trimming steps were applied")

    writer = seq_io.fasta_writer if (args.out_format == 'fasta') else \
        seq_io.fastq_writer

    paired = True if (args.interleaved or args.r_file) else False
 
    if paired:
        print("\nProcessing input as paired-end reads", file=sys.stderr)
        seq_io.logger(args.log, "Record\tForward length\tForward trimmed "
            "length\tReverse length\tReverse trimmed length\n")

        out_s = args.out_s if args.out_s else None
        out_r = out_f if ((args.interleaved or args.out_interleaved) and not \
            args.out_r) else args.out_r

        pairs_passed = discarded_pairs = fsingles = rsingles = 0
        for i, (forward, reverse) in enumerate(iterator):
            identifier = forward['identifier']
            forig = len(forward['sequence'])
            rorig = len(reverse['sequence'])

            forward = apply_trimming(forward, trim_steps, args.qual_type, 
                fheadcrop, fcrop, args.trunc_n)
            ftrim = len(forward['sequence'])

            reverse = apply_trimming(reverse, trim_steps, args.qual_type, 
                rheadcrop, rcrop, args.trunc_n)
            rtrim = len(reverse['sequence'])

            # both good
            if ftrim >= fminlen and rtrim >= rminlen:
                pairs_passed += 1
                writer(out_f, forward)
                writer(out_r, reverse)
            # forward orphaned, reverse filtered
            elif ftrim >= fminlen and rtrim < rminlen:
                fsingles += 1
                writer(out_s, forward)
            # reverse orphaned, forward filtered
            elif ftrim < fminlen and rtrim >= rminlen:
                rsingles += 1
                writer(out_s, reverse)
            # both discarded
            else:
                discarded_pairs += 1

            seq_io.logger(args.log, "{}\t{}\t{}\t{}\t{}\n".format(identifier, 
                forig, ftrim, rorig, rtrim))

        try:
            i += 1
        except UnboundLocalError:
            seq_io.print_error("error: no sequences were found to process")

        total = i * 2
        passed = pairs_passed * 2 + fsingles + rsingles
        print("\nRecords processed:\t{!s} ({!s} pairs)\nPassed filtering:\t"
            "{!s} ({:.2%})\n  Paired reads kept:\t{!s} ({:.2%})\n  Forward "
            "only kept:\t{!s} ({:.2%})\n  Reverse only kept:\t{!s} ({:.2%})"
            "\nRead pairs discarded:\t{!s} ({:.2%})\n".format(total, i,
            passed, passed / total, pairs_passed, pairs_passed / i,
            fsingles, fsingles / total, rsingles, rsingles / total,
            discarded_pairs, discarded_pairs / i), file=sys.stderr)

    else:
        print("\nProcessing input as single-end reads", file=sys.stderr)
        seq_io.logger(args.log, "Record\tLength\tTrimmed length\n")

        if args.out_s:
            print("\nwarning: argument --singles used with single-end reads"
                "... ignoring\n", file=sys.stderr)

        discarded = 0
        for i, record in enumerate(iterator):
            if i == 0:
                first_read = record['identifier']
            elif i == 1:
                if first_read == record['identifier'] and not args.force:
                    seq_io.print_error("warning: the input fastq appears to "
                        "contain interleaved paired-end reads. Please run with "
                        "the --force flag to proceed with processing the data "
                        "as single-end reads")

            origlen = len(record['sequence'])
            record = apply_trimming(record, trim_steps, args.qual_type,
                fheadcrop, fcrop, args.trunc_n)
            trimlen = len(record['sequence'])
            
            if trimlen >= fminlen:
                writer(out_f, record)
            else:
                discarded += 1

            seq_io.logger(args.log, "{}\t{}\t{}\n".format(record['identifier'],
                origlen, trimlen))

        try:
            i += 1
        except UnboundLocalError:
            seq_io.print_error("error: no sequences were found to process. Is "
                "the input properly formatted?")
 
        passed = i - discarded
        print("\nRecords processed:\t{!s}\nPassed filtering:\t{!s} "
        "({:.2%})\nRecords discarded:\t{!s} ({:.2%})\n".format(i, passed,
        passed / i, discarded, discarded / i), file=sys.stderr)

if __name__ == "__main__":
    main()
    sys.exit(0)
