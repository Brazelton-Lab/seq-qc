def fastq_parser(filehandle, line=None):
    """
    Iterator over the given FASTQ file handle returning records (inspired by \
    screeds fastq_iter).
    """
    if line is None:
        line = filehandle.readline()
    while line:
        data = {}
        if not line.startswith('@'):
            raise IOError("Bad FASTQ format: no '@' at beginning of line")

        try:
            data['identifier'], data['description'] = line[1:].strip().split(' ', 1)
        except ValueError:  # No optional annotations
            data['identifier'] = line[1:].strip()
            data['description'] = ''
            pass

        # Extract the sequence lines
        sequence = []
        line = filehandle.readline()
        while not line.startswith('+') and not line.startswith('#'):
            sequence.append(line.strip())
            line = filehandle.readline()

        data['sequence'] = ''.join(sequence)

        # Extract the quality lines
        quality = []
        line = filehandle.readline()
        seqlen = len(data['sequence'])
        aclen = 0
        while not line == '' and aclen <= seqlen:
            quality.append(line.strip())
            aclen += len(line)
            line = filehandle.readline()

        data['quality'] = ''.join(quality)
        if len(data['sequence']) != len(data['quality']):
            raise IOError('sequence and quality strings must be '
                          'of equal length')

        yield data


def fasta_iter(handle, header=None):
    """Iterate over FASTA file and return FASTA entries

    :param handle: FASTA file handle
    :type handle: File Object

    :param header: Header line of entry to start iterating from
    :type header: str

    :return: dictionary containing FASTA data
    :rtype: dict

    @Alex Hyer
    """

    strip = str.strip  # Speed trick: reduces function calls

    if header is None:
        header = handle.next()  # Read first FASTA entry header
    else:
        header = header  # Set header to given header

    try:  # Manually construct a for loop to improve speed by using 'next'

        while True:  # Loop until StopIteration Exception raised

            line = handle.next()

            data = []

            if not header[0] == '>':
                raise IOError('Bad FASTA format: no ">" at beginning of line')

            header = strip(header)
            try:
                data['identifier'], data['description'] = strip(header)[1:].split(' ', 1)
            except ValueError:  # No description
                data['identifier'] = strip(header)[1:]
                data['description'] = ''

            # Collect sequence lines into a list
            sequence_list = []
            while line and not line[0] == '>':
                sequence_list.append(strip(line))
                line = handle.next()  # Raises StopIteration on last line
            header = line  # Store current line so it's not lost next iteration

            data['sequence'] = ''.join(sequence_list)
            yield data

    except StopIteration:
        pass
    finally:  # Yield last FASTA entry
        data['sequence'] = ''.join(sequence_list)
        yield data


def interleaved_parser(file_iter):
    """Read pairs from a stream (inspired by khmer).

    A generator that yields singletons pairs from a stream of FASTA/FASTQ 
    records.  Yields (n, is_pair, r1, r2) where 'r2' is None if is_pair is 
    False.

    Usage::

       for read1, read2 in interleaved_iter(...):
          ...
    """
    import pairs

    record = None
    prev_record = None

    # handle the majority of the stream.
    for record in file_iter:
        if prev_record:
            if pairs.verify_paired(prev_record, record):
                yield (prev_record, record) # it's a pair!
                record = None
            else: # orphan.
                raise pairs.UnpairedReadsError("Unpaired reads found. Data "
                    "may contain orphans or is not ordered properly", 
                    prev_record, record)

        prev_record = record
        record = None
