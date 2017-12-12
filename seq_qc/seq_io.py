from __future__ import print_function

from arandomness.argparse import Open
from bio_utils.iterators import fasta_iter, fastq_iter
import bz2
import gzip
import io
from seq_qc import pairs
import sys
import textwrap


class Paired:
    """A wrapper to FastaEntry and FastqEntry classes for paired reads. 

    Attributes:
        forward (class): first read of the pair
        
        reverse (class): second read of the pair
    """
    def __init__(self, forward, reverse):
        """Attributes to store paired entry data"""

        self.forward = forward
        self.reverse = reverse


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
        return handle_interleaved(f_iter)
    elif reverse:
        r_iter = parser(reverse)
        return handle_split(izip(f_iter, r_iter))
    else:
        return (f_iter)


def handle_interleaved(file_iter):
    """Read interleaved pairs from a stream (inspired by khmer).

    A generator that yields singletons pairs from a stream of FASTA/FASTQ
    records.  Yields (r1, r2) where 'r2' is None if is_pair is
    False.

    Usage::

       for read1, read2 in interleaved_wrapper(...):
          ...
    """

    record = None
    prev_record = None

    # Handle the majority of the stream.
    for record in file_iter:
        if prev_record:
            if pairs.verify_paired(prev_record, record):
                yield Paired(prev_record, record)  #records are paired
                record = None
            else:  #one of the pairs is orphaned
                raise pairs.UnpairedReadsError("Unpaired reads found. Data "
                    "may contain orphans or is not ordered properly",
                    prev_record, record)

        prev_record = record
        record = None


def handle_split(file_iter):
    """Read split pairs from a stream
    """
    for record1, record2 in file_iter:
        if pairs.verify_paired(record1, record2):
            yield Paired(record1, record2)  #records are paired
        else:
            raise pairs.UnpairedReadsError("Unpaired reads found. Verify that "
                "the order of the forward and reverse reads is the same", 
                record1, record2)


def print_error(message):
    print(textwrap.fill(message, 79), file=sys.stderr)
    sys.exit(1)


def program_info(prog, args, version):
    print("{} {!s}".format(prog, version), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}".format(' '.join(args)), 79), file=sys.stderr)
