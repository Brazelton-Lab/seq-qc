from __future__ import print_function

from arandomness.argparse import Open
from bio_utils.iterators import fasta_iter, fastq_iter
import bz2
import gzip
import io
import sys
import textwrap


def read_iterator(forward, reverse=None, interleaved=False, f_format='fastq'):
    """Generates an object to iterate over. If reads are paired-end, each 
    record will contain both the forward and reverse sequences of the pair.
    """
    try:
        # Python2
        from itertools import izip
    except ImportError:
        # Python3
        izip = zip

    if f_format == 'fasta':
        parser = fasta_iter
    elif f_format == 'fastq':
        parser = fastq_iter

    f_iter = parser(forward)

    # Wrap pairs if required
    if interleaved:
        return interleaved_wrapper(f_iter)
    elif reverse:
        r_iter = parser(reverse)
        return izip(f_iter, r_iter)
    else:
        return f_iter


def interleaved_wrapper(file_iter):
    """Read pairs from a stream (inspired by khmer).

    A generator that yields singletons pairs from a stream of FASTA/FASTQ
    records.  Yields (r1, r2) where 'r2' is None if is_pair is
    False.

    Usage::

       for read1, read2 in interleaved_wrapper(...):
          ...
    """
    from seq_qc import pairs

    record = None
    prev_record = None

    # Handle the majority of the stream.
    for record in file_iter:
        if prev_record:
            if pairs.verify_paired(prev_record, record):
                yield (prev_record, record) #it's a pair!
                record = None
            else: # orphan.
                raise pairs.UnpairedReadsError("Unpaired reads found. Data "
                    "may contain orphans or is not ordered properly",
                    prev_record, record)

        prev_record = record
        record = None


def print_error(message):
    print(textwrap.fill(message, 79), file=sys.stderr)
    sys.exit(1)


def program_info(prog, args, version):
    print("{} {!s}".format(prog, version), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}".format(' '.join(args)), 79), file=sys.stderr)
