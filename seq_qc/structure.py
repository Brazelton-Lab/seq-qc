from seq_qc.seq_io import print_error

def reverse_complement(seq):
    """Take the reverse complement of a nucleotide sequence"""
    r = reverse(seq)
    rc = complement(r)

    return rc


__complementaryBasePairs = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

def complement(seq):
    """Take the complement of a nucleotide sequence"""
    try:
        c = "".join(__complementaryBasePairs[i] for i in seq)
    except KeyError:
        print_error("error: non-cannonical representation of a nucleotide "
                    "sequence provided. String is {0}".format(seq))

    return c


def reverse(seq):
    """Take the reverse of a sequence"""
    return seq[::-1]


def reverse_complement_paired(seq1, seq2):
    """
    Return a tuple containing the reverse-complement of paired reads
    """
    r = reverse_paired(seq1, seq2)
    return tuple([complement(i) for i in r])


def reverse_paired(seq1, seq2):
    """
    Return a tuple containing the reverse of the paired reads
    """
    return (seq2, seq1)
