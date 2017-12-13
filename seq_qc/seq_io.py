from __future__ import print_function

import argparse
from bio_utils.iterators import fasta_iter, fastq_iter
from seq_qc import pairs
import sys
import textwrap


class Open(argparse.Action):
    """Argparse Action that detects and opens compressed files for rw

    Attributes:
        option_strings (list): list of str giving command line flags that
                               call this action

        dest (str): Namespace reference to value

        mode (str): mode to pass to (de)compression algorithm

        nargs (bool): True if multiple arguments specified

        **kwargs (various): optional arguments to pass to argparse and algo
    """

    def __init__(self, option_strings, dest, mode='rb', nargs=None, **kwargs):
        """Initialize class and spawn self as Base Class w/o nargs

        Warns:
            ImportError: if Open cannot import a compression library,
                         it warns the user that it cannot open the
                         corresponding file type

        Raises:
            ValueError: if nargs is not None, Open does not accept nargs
        """
        from importlib import import_module
        from warnings import warn
        from os import linesep

        # Only accept a single value to analyze
        if nargs is not None:
            raise ValueError('nargs not allowed for Open')

        # Call self again but without nargs
        super(Open, self).__init__(option_strings, dest, **kwargs)

        # Store and establish variables used in __call__
        self.kwargs = kwargs
        self.mode = mode.lower().strip()
        self.modules = {}

        modules_to_import = {
            'bz2': 'BZ2File',
            'gzip': 'GzipFile',
            'lzma': 'LZMAFile'
        }

        # Dynamically import compression libraries and warn about failures
        for mod, _class in modules_to_import.items():
            try:
                self.modules[_class] = getattr(import_module(mod), _class)
            except (ImportError, AttributeError) as e:
                self.modules[_class] = open
                warn('Cannot process {0} files due to following error:'
                     '{1}{2}{1}You will need to install the {0} library to '
                     'properly use these files. Currently, such files will '
                     'open in text mode.'.format(mod, linesep, e))
    # Credits: https://stackoverflow.com/questions/13044562/
    # python-mechanism-to-identify-compressed-file-type-and-uncompress
    def __call__(self, parser, namespace, value, option_string=None, **kwargs):
        """Detects and opens compressed files

        Args:
            parser (ArgumentParser): parser used to generate values

            namespace (Namespace): namespace to set values for

            value (str): actual value specified by user

            option_string (str): argument flag used to call this function

            **kwargs (various): optional arguments later passed to the
                                compression algorithm
        """
        from inspect import getfullargspec
        import io

        filename = value  # For readability

        algo = io.open  # Default to plaintext

        # Capture any mode that isn't read, such as write or append
        if self.mode.lstrip('U')[0] != 'r':

            algo_map = {
                'bz2': self.modules['BZ2File'],
                'gz':  self.modules['GzipFile'],
                'xz':  self.modules['LZMAFile']
            }

            # Base compression algorithm on file extension
            ext = value.split('.')[-1]
            try:
                algo = algo_map[ext]
            except KeyError:
                pass

        # Basically read mode
        else:

            file_sigs = {
                b'\x42\x5a\x68': self.modules['BZ2File'],
                b'\x1f\x8b\x08': self.modules['GzipFile'],
                b'\xfd7zXZ\x00': self.modules['LZMAFile']
                }

            max_len = max(len(x) for x in file_sigs.keys())


            # Check if file is connected to a tty device
            if sys.stdin.isatty():
                # Check beginning of file for signature
                with io.open(filename, 'rb') as in_handle:
                    start = in_handle.read(max_len)
                    for sig in file_sigs.keys():
                        if start.startswith(sig):
                            algo = file_sigs[sig]
                            break


        # Filter all **kwargs by the args accepted by the compression algo
        algo_args = set(getfullargspec(algo).args)
        good_args = set(self.kwargs.keys()).intersection(algo_args)
        _kwargs = {arg: self.kwargs[arg] for arg in good_args}


        # Open the file using parameters defined above and store in namespace
        try:
            handle = algo(value, mode=self.mode, **_kwargs)
        except ValueError:
            mode = self.mode.lstrip('U')[0]
            handle = io.TextIOWrapper(algo(value, mode=mode, **_kwargs))

        setattr(namespace, self.dest, handle)


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
            raise pairs.UnpairedReadsError("Unpaired reads found. Data may "
                "contain orphans or is not ordered properly",
                record1, record2)


def print_error(message):
    print(textwrap.fill(message, 79), file=sys.stderr)
    print("", file=sys.stderr)
    sys.exit(1)


def print_warning(message):
    print(textwrap.fill(message, 79), file=sys.stderr)
    print("", file=sys.stderr)


def program_info(prog, args, version):
    print("{} {!s}".format(prog, version), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}".format(' '.join(args)), 79), file=sys.stderr)
    print("", file=sys.stderr)
