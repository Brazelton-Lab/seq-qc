from __future__ import print_function

import bz2
import gzip
import io
from itertools import izip
import seq_parsers
import sys
import textwrap

def open_output(filename):
    """
    Decide how to open an output file for writing based on the file extension
    """
    extension = filename.split('.')[-1]
    if extension == 'gz':
        return gzip.GzipFile(filename, 'wb')
    elif extension == 'bz2':
        return bz2.BZ2File(filename, 'wb')
    else:
        return open(filename, 'w')

def open_input(filename):
    """
    Make a best-effort guess as to how to open/parse the given sequence file.

    Deals with .gz, FASTA, and FASTQ records.
    """
    file_signatures = {
        "\x1f\x8b\x08": "gz",
        "\x42\x5a\x68": "bz2",
        # "\x50\x4b\x03\x04": "zip"
    }  # Inspired by http://stackoverflow.com/a/13044946/1585509
    try:
        bufferedfile = io.open(file=filename, mode='rb', buffering=8192)
    except IOError:
        print_message("error: unable to open {} for reading".format(filename),
            sys.stderr)
    num_bytes_to_peek = max(len(x) for x in file_signatures)
    file_start = bufferedfile.peek(num_bytes_to_peek)
    compression = None
    for signature, ftype in file_signatures.items():
        if file_start.startswith(signature):
            compression = ftype
            break

    if compression is 'bz2':
        filehandle = bz2.BZ2File(filename)
        peek = filehandle.peek(1)
    elif compression is 'gz':
        if not bufferedfile.seekable():
            print_message("gziped data not streamable, pipe through zcat \
                             first", sys.stderr)
        peek = gzip.GzipFile(filename).read(1)
        filehandle = gzip.GzipFile(filename)
    else:
        peek = bufferedfile.peek(1)
        filehandle = bufferedfile

    parser = None
    try:
        if peek[0] == '>':
            parser = seq_parsers.fasta_parser
        elif peek[0] == '@':
            parser = seq_parsers.fastq_parser
    except IndexError as err:
        return []  # empty file

    if parser is None:
        # return the filehandle instead of an iterator
        return filehandle
    else:
        return parser(filehandle)

def get_iterator(forward, reverse=None, interleaved=False):
    f_iter = open_input(forward)
    if interleaved:
        iterator = seq_parsers.interleaved_parser(f_iter)
    elif reverse:
        r_iter = open_input(reverse)
        iterator = izip(f_iter, r_iter)
    else:
        iterator = f_iter

    return iterator

def print_message(message, destination=sys.stdout):
    print(textwrap.fill(message, 79), file=destination)
    if destination == sys.stderr:
        sys.exit(1)

def start_message(prog, args, version):
    print("{} {!s}".format(prog, version))
    print_message("Starting {} with arguments: {}".format(prog, ' '.join(args)))

def logger(loghandle, output):
    if loghandle:
        loghandle.write(output)

def fastq_writer(filehandle, record):
    """
    Write a record to a fastq file
    """
    try:
        quality = record['quality']
    except KeyError:
        print_message("error: could not find quality score information for "
            "the sequences. Please verfiy that the input file is in fastq "
            "format.", sys.stderr)

    output = ("@{} {}\n{}\n+\n{}\n".format(record['identifier'],
        record['description'], record['sequence'], quality))
    if filehandle:
        filehandle.write(output)

def fasta_writer(filehandle, record):
    """
    Write a record to a fasta file
    """
    output = (">{} {}\n{}\n".format(record['identifier'], 
        record['description'], record['sequence']))
    if filehandle:
        filehandle.write(output)
