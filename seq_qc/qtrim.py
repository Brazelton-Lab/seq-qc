#! /usr/bin/env python
"""
Standard sequence quality control tool that can be used for base cropping, 
trimming sequences by quality score, and filtering reads that fail to meet
a minimum length threshold post cropping and trimming.

For single-end and interleaved reads:
    qtrim [options] [-o output] input
 
For split paired-end reads:
    qtrim [option] -o out.forward -v out.reverse -s out.singles -r in.reverse
        in.forward

Input files must be in FASTQ format. Output can be in either FASTQ or FASTA 
format. Compression using gzip and bzip2 algorithms is automatically detected 
for input files. To compress output, add the appropriate file extension to the 
file names (.gz, .bz2). For single-end or interleaved reads, use /dev/stdin to 
indicate that input should be taken from standard input (stdin). Similarly, 
leaving out the -o argument will cause output to be sent to standard output 
(stdout).
"""

from __future__ import division
from __future__ import print_function

from arandomness.argparse import CheckThreads, Open
import argparse
from bio_utils.iterators import FastqEntry
from multiprocessing import cpu_count, Lock, Process, Queue, Value
from seq_qc.pairs import verify_paired
from seq_qc import seq_io, trim
from subprocess import check_output
import sys
from time import sleep, time

__author__ = "Christopher Thornton"
__license__ = 'GPLv2'
__maintainer__ = 'Christopher Thornton'
__version__ = "2.0.0"
__status__ = "Production"


class Counter(object):
    """A synchronized shared counter to use with the multiprocessing module.

    Credit: Eli Benderksy
    """
    def __init__(self, initval=0):
        self.val = Value('i', initval)
        self.lock = Lock()

    def increment(self):
        with self.lock:
            self.val.value += 1

    def value(self):
        with self.lock:
            return self.val.value


def trim_reads(rqueue, wqueue, steps, trunc_n, crop, hcrop, minlen, p, d, s1, s2, offset=33):
    """Processes accepting input from the read queue are responsible for 
    performing trimming based on quality score. Trimmed reads will be put into 
    the write queue for subsequent threshold evaluation and writing.

    Args:
         rqueue (Queue): multiprocessing Queue class containing records to 
                         process

         wqueue (Queue): multiprocessing Queue class to input trimmed records 
                         designated for writing

         steps (dict): dictionary of trimming steps to apply, with trimming 
                       functions as keys and user-supplied arguments as the 
                       values

         trunc_n (function): function for trimming reads to position of first 
                             ambiguous base, or returning data unchanged

         crop (array): array storing indices to crop sequences to

         hcrop (array): array storing indices where sequences should start 
                        from

         offset (int): quality score offset. Can be 33 (Sanger) or 64
    """

    # Loop until queue contains kill message
    while True:

        entry = rqueue.get()

        # Break on kill message
        if entry == 'DONE':
            break

        try:
            entry = (entry.forward, entry.reverse)
        except AttributeError:
            records = (entry,)

        trimmed = []
        for j, record in enumerate(records):
            scores = trim.translate_quality(record.quality, offset)

            seqlen = len(scores)
            lastbase = seqlen if not crop[j] or (crop[j] > seqlen) else crop[j]

            start = hcrop[j]  #initial starting index of a sequence
            end = 0  #intial number of bases to remove from sequence end

            for step, value in steps:
                newstart, newend = step(scores[start: lastbase - end], value)
                start += newstart
                end += newend

            last = seqlen - end

            record.sequence, record.quality = record.sequence[start: last], \
                record.quality[start: last]

            trimmed.append(trunc_n(record))

        trimlen = [len(i.sequence) for i in trimmed]

        first_greater = trimlen[0] >= minlen[0]
        try:
            second_greater = trimlen[1] >= minlen[1]

        except IndexError:
            # Record passed length threshold
            if first_greater:
                p.increment()
                wqueue.put((("forward", trimmed[0]),))

            else:
                d.increment()

        else:
            # Both read pairs passed length threshold
            if first_greater and second_greater:
                p.increment()
                wqueue.put((("forward", trimmed[0]), ("reverse", trimmed[1])))

            # Forward reads orphaned, reverse reads failed length threshold
            elif trimlen[0] >= minlen[0] and trimlen[1] < minlen[1]:
                s1.increment()
                wqueue.put((("forward", trimmed[0]),))

            # Reverse reads orphaned, forward reads failed length threshold
            elif trimlen[0] < minlen[0] and trimlen[1] >= minlen[1]:
                s2.increment()
                wqueue.put((("reverse", trimmed[1]),))

            # Both read pairs failed length threshold and were discarded
            else:
                d.increment()

    # Send kill message to thread responsible for writing
    wqueue.put('DONE')


def write_reads(queue, fwrite, rwrite, swrite):
    """Processes accepting input from the write queue are responsible for 
    evaluating the results of the length comparison and for writing records
to the output stream.

    Args:
         queue (Queue): multiprocessing Queue class containing trimmed records

         fwrite (function): function for writing records to the forward file
         
         rwrite (function): function for writing records to the reverse file

         swrite (function): function for writing records to the singles file
    """

    mapper = {"singles": swrite, "forward": fwrite, "reverse": rwrite}

    while True:

        trimmed = queue.get()

        # Break on kill message
        if trimmed == 'DONE':
            break
        for record in trimmed:
            print(record[1].write(), end='', file=mapper[record[0]])


def parse_colons(argument):
    try:
        window, score = argument.split(':')
    except ValueError:
        seq_io.print_error("error: the input provided to sliding-window is "
                           "formatted incorrectly. See --help for usage")
    else:
        if score.isdigit():
            score = int(score)
        else:
            seq_io.print_error("error: the quality score threshold provided "
                               "to sliding-window must be an integer value")
        if window.isdigit():
            window = int(window)
        else:
            try:
                window = float(window)
            except ValueError:
                seq_io.print_error("error: the window-size provided to "
                                   "sliding-window must be either an integer "
                                   "value or a fraction")

    return (window, score)


def parse_commas(args, argname):
    args = [i.lstrip() for i in args.split(",")]

    if 1> len(args) > 2:
        seq_io.print_error("error: only one or two integer values should be "
            "provided to {0}".format(argname))

    try:
        arg1 = int(args[0])
        arg2 = int(args[1])
    except ValueError:
        seq_io.print_error("error: input to {0} must be one or more integer "
                           "values in the form INT or INT,INT".format(argname))
    except IndexError:
        arg1 = arg2 = int(args[0])

    return (arg1, arg2)


def do_nothing(*args):
    pass


def as_is(args):
    return args


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fhandle', 
        metavar='in1.fastq',
        type=str,
        action=Open,
        mode='rb',
        default=sys.stdin,
        help="input reads in fastq format. Can be a file containing either "
             "single-end or forward/interleaved reads if reads are paired-end "
             "[required]")
    input_arg = parser.add_mutually_exclusive_group(required=False)
    input_arg.add_argument('--interleaved',
        action='store_true',
        help="input is interleaved paired-end reads")
    input_arg.add_argument('-r', '--reverse',
        dest='rhandle', 
        metavar='in2.fastq', 
        action=Open,
        mode='rb',
        help="input reverse reads in fastq format")
    parser.add_argument('-o', '--out', 
        metavar='FILE', 
        dest='out_f',
        type=str,
        action=Open,
        mode='w',
        default=sys.stdout,
        help="output trimmed reads [default: stdout]")
    output_arg = parser.add_mutually_exclusive_group(required=False)
    output_arg.add_argument('-v', '--out-reverse', 
        metavar='FILE', 
        dest='out_r',
        type=str,
        action=Open,
        mode='w',
        help="output trimmed reverse reads")
    output_arg.add_argument('--out-interleaved', 
        dest='out_interleaved',
        action='store_true',
        help="output interleaved paired-end reads, even if input is split")
    parser.add_argument('-s', '--singles', 
        metavar='FILE', 
        dest='out_s',
        type=str,
        action=Open,
        mode='w',
        help="output trimmed orphaned reads")
    parser.add_argument('-q', '--qual-offset', 
        metavar='TYPE', 
        dest='offset',
        type=int, 
        choices=[33, 64],
        default=33,
        help="ASCII base quality score encoding [default: 33]. Options are "
             "33 (for phred33) or 64 (for phred64)")
    parser.add_argument('-m', '--min-len', 
        metavar='LEN [,LEN]', 
        dest='minlen',
        type=str, 
        help="filter reads shorter than the minimum length threshold [default:"
             " 0]. Different values can be provided for the forward and "
             "reverse reads, respectively, by separating them with a comma "
             "(e.g. 80,60), or a single value can be provided for both")
    trim_args = parser.add_argument_group('trimming options')
    trim_args.add_argument('-O', '--trim-order', 
        metavar='ORDER',
        dest='trim_order',
        type=str,
        default='ltw',
        help="order in which the trimming methods should be applied. "
             "Available methods are l (leading), t (trailing), and w "
             "(sliding-window) [default: ltw]")
    trim_args.add_argument('-W', '--sliding-window', 
        metavar='FRAME',
        dest='sw',
        type=parse_colons,
        help="trim read ends using a sliding window approach. Input should be "
             "of the form 'window_size:qual_threshold', where 'qual_threshold' "
             "is an integer between 0 and 42 and 'window_size' can either be "
             "length in bases or fraction of total read length")
    trim_args.add_argument('-H', '--headcrop', 
        metavar='INT [,INT]',
        type=str,
        help="remove exactly the number of bases specified from the start of "
             "the reads [default: 0]. Different values can be provided for "
             "the forward and reverse reads, respectively, by separating them "
             "with a comma (e.g. 2,0), or a single value can be provided for "
             "both. Cropping will always be applied first")
    trim_args.add_argument('-C', '--crop', 
        metavar='INT [,INT]',
        type=str,
        help="crop reads to the specified position [default: off]. The "
             "value(s) should be less than the maximum read length, otherwise "
             "no cropping will be applied. Different values can be provided "
             "for the forward and reverse reads, respectively, by separating "
             "them with a comma (e.g. 120,115), or a single value can be "
             "provided for both. Cropping will always be applied first")
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
        help="truncate sequence at position of first ambiguous base [default: "
             "off]. Truncation will always be applied last.")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    parser.add_argument('-t', '--threads',
        action=CheckThreads,
        type=int,
        default=1,
        help='number of threads to use for trimming [default: 1]')
    args = parser.parse_args()
    all_args = sys.argv[1:]

    seq_io.program_info('qtrim', all_args, __version__)


    # Fail if insufficient directive supplied
    if args.rhandle and not (args.out_r or args.out_interleaved):
        parser.error("one of -v/--out-reverse or --out-interleaved is required "
            "when the argument -r/--reverse is used")


    # Track program run-time
    start_time = time()


    # Assign variables based on arguments supplied by the user
    crop = parse_commas(args.crop, "crop") if args.crop else (None, None)
    hcrop = parse_commas(args.headcrop, "headcrop") if \
        args.headcrop else (0, 0)
    minlen = parse_commas(args.minlen, "minlen") if args.minlen \
        else (0, 0)
    out_f = args.out_f
    paired = True if (args.interleaved or args.rhandle) else False
    trunc_n = trim.truncate_by_n if args.trunc_n else as_is


    # Prepare the iterator based on dataset type
    iterator = seq_io.read_iterator(args.fhandle, args.rhandle, \
                                    args.interleaved, "fastq")


    # Populate list of trimming tasks to perform on reads
    trim_tasks = {'l': (trim.trim_leading, args.lead_score), 
        't': (trim.trim_trailing, args.trail_score), 
        'w': (trim.adaptive_trim, args.sw)}

    trim_steps = []
    for task in args.trim_order:
        value = trim_tasks[task][-1]
        if value:
            trim_steps.append(trim_tasks[task])
    if len(trim_steps) < 1 and not (args.crop or args.headcrop):
        seq_io.print_error("error: no trimming steps were specified")


    # Counters
    discarded = Counter(0)
    passed = Counter(0)

    # Check dataset type (paired or single-end) 
    if paired:
        print("\nProcessing input as paired-end reads", file=sys.stderr)

        out_s = args.out_s if args.out_s else do_nothing
        out_r = out_f if ((args.interleaved or args.out_interleaved) and not \
            args.out_r) else args.out_r

        output = "\nRecords processed:\t{!s}\nPassed filtering:\t{!s} " \
                 "({:.2%})\n  Reads pairs kept:\t{!s} ({:.2%})\n  Forward " \
                 "only kept:\t{!s} ({:.2%})\n  Reverse only kept:\t{!s} " \
                 "({:.2%})\nRecords discarded:\t{!s} ({:.2%})\n"

        singles1 = Counter(0)
        singles2 = Counter(0)

    else:
        print("\nProcessing input as single-end reads", file=sys.stderr)

        if args.out_s:
            print("\nwarning: argument --singles used with single-end reads"
                  "... ignoring\n", file=sys.stderr)

        out_s = do_nothing
        out_r = do_nothing

        output = "\nRecords processed:\t{!s}\nPassed filtering:\t{!s} ({:.2%})" \
                 "\nRecords discarded:\t{!s} ({:.2%})\n"

        singles1 = singles2 = None

    max_read_threads = args.threads - 1 if args.threads > 1 else 1
    read_queue = Queue(max_read_threads)  # Max queue prevents race conditions
    write_queue = Queue()

   
    # Initialize threads to process reads and writes
    read_processes = []
    for i in range(max_read_threads):
        read_processes.append(Process(target=trim_reads, args=(read_queue, write_queue, trim_steps, trunc_n, crop, hcrop, minlen, passed, discarded, singles1, singles2, args.offset,)))
        read_processes[i].start()

    write_process = Process(target=write_reads, args=(write_queue, out_f, out_r, out_s,))
    write_process.start()


    # Iterate over reads, populating read queue for trimming
    for i, records in enumerate(iterator):
        read_queue.put(records)


    # Send a kill message to each thread designated for reading via read_queue
    for j in read_processes:
        read_queue.put('DONE')


    # Wait for processes to finish before continuing
    for process in read_processes:
        process.join()

    write_process.join()


    # Make sure input file non-empty
    try:
        i += 1
    except UnboundLocalError:
        seq_io.print_error("error: no sequences were found to process")


    # Calculate and print output statistics
    p = passed.value()
    d = discarded.value()

    if paired:
        i = i * 2
        s1, s2 = singles1.value(), singles2.value()
        pairs = p
        total = (p * 2 + s1 + s2)
        frac_pairs = (pairs * 2) / i
        frac_s1 = s1 / i
        frac_s2 = s2 /i
        d = (i - d * 2) - (s1 + s2)
    else:
        s1 = s2 = frac_s1 = frac_s2 = pairs = frac_pairs = None
        total = p

    frac_d = d / i
    frac_total = total / i
    stats = [i, total, frac_total] + [i for i in (pairs, frac_pairs, s1, \
             frac_s1, s2, frac_s2) if i != None] + [d, frac_d]
    print(output.format(*tuple(stats)), file=sys.stderr)


    # Calculate and print program run-time
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("It took {:.2e} minutes to process {!s} records"\
          .format(total_time, i), file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
