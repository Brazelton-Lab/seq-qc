#! /usr/bin/env python
"""
Sequence quality control tool for common preprocessing tasks such as cropping, 
quality-based trimming, and minimum length filtering.

For single-end and interleaved reads:
    qtrim [options] [-o out.reads] in.reads
 
For split paired-end reads:
    qtrim [option] [-o out.forward] [-v out.reverse] [-s out.singles] 
        -r in.reverse in.forward

FASTQ is the only supported file format. The compression algorithm is 
automatically detected for input files through their file extension. To 
compress the output files, add the appropriate file extension to the file names 
(.gz, .bz2). For single-end or interleaved reads, use /dev/stdin to indicate 
that input is from standard input (stdin). Similarly, leaving out '--output' 
will result in the output being sent to standard output (stdout).
"""

from __future__ import division
from __future__ import print_function

from arandomness.argparse import CheckThreads
import argparse
from multiprocessing import cpu_count, Lock, Process, Queue, Value
from seq_qc import seq_io, trim
from subprocess import check_output
import sys
from time import sleep, time

__author__ = "Christopher Thornton"
__license__ = 'GPLv2'
__maintainer__ = 'Christopher Thornton'
__status__ = "Production"
__version__ = "2.0.1"


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


def trim_reads(rqueue, wqueue, steps, trunc_n, crop, hcrop, offset=33):
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
    # Loop until read queue contains kill message
    while True:

        entry = rqueue.get()

        # Break on kill message
        if entry == 'DONE':
            break

        try:
            records = (entry.forward, entry.reverse)
        except AttributeError:
            records = (entry,)

        trimmed = []
        for i, record in enumerate(records):
            scores = trim.translate_quality(record.quality, offset)

            seqlen = len(scores)
            lastbase = seqlen if not crop[i] or (crop[i] > seqlen) else crop[i]

            start = hcrop[i]  #initial starting index of a sequence
            end = 0  #intial number of bases to remove from sequence end

            for step, value in steps:
                newstart, newend = step(scores[start: lastbase - end], value)
                start += newstart
                end += newend

            last = lastbase - end

            record.sequence, record.quality = record.sequence[start: last], \
                record.quality[start: last]

            trimmed.append(trunc_n(record))

        wqueue.put(trimmed)


def write_reads(queue, fout, rout, sout, minlen, p, d, s1, s2):
    """Processes accepting input from the write queue are responsible for 
    evaluating the results of the length comparison and for writing records
    to the output stream.

    Args:
         queue (Queue): multiprocessing Queue class containing trimmed records

         fwrite (function): function for writing records to the forward file
         
         rwrite (function): function for writing records to the reverse file

         swrite (function): function for writing records to the singles file

         minlen(int): threshold for determining whether to write a record
         
         p, d, s1, s2 (Counter): counter class for synchronizing incremention
                                 between processes
    """
    lock = Lock()

    fwrite = fout.write
    rwrite = rout.write if rout else do_nothing
    swrite = sout.write if sout else do_nothing

    # Loop until write queue contains kill message
    while True:

        trimmed = queue.get()

        # Break on kill message
        if trimmed == 'DONE':
            break

        try:
            trimlen = [len(i.sequence) for i in trimmed]
        except AttributeError:
            print_error("error: ")

        first_greater = (trimlen[0] >= minlen[0])
        try:
            second_greater = (trimlen[1] >= minlen[1])
        except IndexError:
            # Record passed length threshold
            if first_greater:
                p.increment()
                fwrite(trimmed[0].write())
            else:
                d.increment()
        else:
            # Both read pairs passed length threshold
            if first_greater and second_greater:
                p.increment()
                with lock:
                    fwrite(trimmed[0].write())
                    rwrite(trimmed[1].write())

            # Forward reads orphaned, reverse reads failed length threshold
            elif first_greater and trimlen[1] < minlen[1]:
                s1.increment()
                swrite(trimmed[0].write())

            # Reverse reads orphaned, forward reads failed length threshold
            elif trimlen[0] < minlen[0] and second_greater:
                s2.increment()
                swrite(trimmed[1].write())

            # Both read pairs failed length threshold and were discarded
            else:
                d.increment()

    for handle in (fout, rout, sout):
        try:
            handle.close()
        except AttributeError:
            continue


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


def do_nothing(args):
    pass


def self(args):
    return args


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('fhandle', 
        metavar='in1.fastq',
        type=str,
        action=seq_io.Open,
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
        action=seq_io.Open,
        mode='rb',
        help="input reverse reads in fastq format")
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
        help="output trimmed reverse reads")
    parser.add_argument('-s', '--singles', 
        metavar='FILE', 
        dest='out_s',
        type=str,
        action=seq_io.Open,
        mode='wt',
        help="output trimmed orphaned reads")
    parser.add_argument('-q', '--qual-offset', 
        metavar='TYPE', 
        dest='offset',
        type=int, 
        choices=[33, 64],
        default=33,
        help="ASCII base quality score encoding [default: 33]. Options are "
             "33 (phred33) or 64 (phred64)")
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
        help="order that the trimming methods should be applied [default: ltw]"
             ". Available methods are l (leading), t (trailing), and w "
             "(sliding-window)")
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
             "off]. Truncation will always be applied last")
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
    trunc_n = trim.truncate_by_n if args.trunc_n else self


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


    # Counters for trimming statistics
    discarded = Counter(0)
    passed = Counter(0)


    # Assign variables based on dataset type (paired or single-end) 
    if paired:
        print("Processing input as paired-end reads", file=sys.stderr)

        out_s = args.out_s if args.out_s else None
        out_r = out_f if not args.out_r else args.out_r

        output = "\nRecords processed:\t{!s}\nPassed filtering:\t{!s} " \
                 "({:.2%})\n  Reads pairs kept:\t{!s} ({:.2%})\n  Forward " \
                 "only kept:\t{!s} ({:.2%})\n  Reverse only kept:\t{!s} " \
                 "({:.2%})\nRecords discarded:\t{!s} ({:.2%})\n"

        singles1 = Counter(0)
        singles2 = Counter(0)

    else:
        if args.out_s:
            print("warning: argument --singles used with single-end reads"
                  "... ignoring\n", file=sys.stderr)

        if args.out_r:
            print("warning: argument --out-reverse used when input is "
                  "single-end... ignoring\n", file=sys.stderr)

        print("Processing input as single-end reads", file=sys.stderr)

        out_s = None
        out_r = None

        output = "\nRecords processed:\t{!s}\nPassed filtering:\t{!s} ({:.2%})" \
                 "\nRecords discarded:\t{!s} ({:.2%})\n"

        singles1 = singles2 = None


    max_read_threads = args.threads - 1 if args.threads > 1 else 1
    read_queue = Queue(max_read_threads)  # Max queue prevents race conditions
    write_queue = Queue()

   
    # Initialize threads to process reads and writes
    read_processes = []
    for i in range(max_read_threads):
        read_processes.append(Process(target=trim_reads, args=(read_queue, \
            write_queue, trim_steps, trunc_n, crop, hcrop, args.offset,)))
        read_processes[i].start()

    write_process = Process(target=write_reads, args=(write_queue, out_f, \
        out_r, out_s, minlen, passed, discarded, singles1, singles2,))
    write_process.start()


    # Iterate over reads, populating read queue for trimming
    for processed_total, records in enumerate(iterator):
        read_queue.put(records)


    # Send kill message to threads responsible for trimming
    for process in read_processes:
        read_queue.put('DONE')


    # Wait for processes to finish before continuing
    for process in read_processes:
        process.join()


    # Send kill message to threads responsible for trimming
    while not write_queue.empty():
        sleep(1)
    write_queue.put('DONE')


    # Wait for write processes to finish before continuing
    write_process.join()


    # Verify input file non-empty
    try:
        processed_total += 1
    except UnboundLocalError:
        seq_io.print_error("error: no sequences were found to process")


    # Calculate and print output statistics
    p = passed.value()
    d = discarded.value()

    if paired:
        processed_total = processed_total * 2
        s1, s2 = singles1.value(), singles2.value()
        passed_total = (p * 2 + s1 + s2)
        discarded_total = d * 2 + s1 + s2
        pairs = p
        frac_pairs = (pairs * 2) / processed_total
        frac_s1 = s1 / processed_total
        frac_s2 = s2 / processed_total
    else:
        passed_total = p
        discarded_total = d
        processed_total = discarded_total + passed_total
        s1 = s2 = frac_s1 = frac_s2 = pairs = frac_pairs = None

    frac_discarded = discarded_total / processed_total
    frac_passed = passed_total / processed_total
    stats = [processed_total, passed_total, frac_passed] + [i for i in (pairs, frac_pairs, s1, \
             frac_s1, s2, frac_s2) if i != None] + [discarded_total, frac_discarded]
    print(output.format(*tuple(stats)), file=sys.stderr)


    # Calculate and print program run-time
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("It took {:.2e} minutes to process {!s} records\n"\
          .format(total_time, processed_total), file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
