from __future__ import division
from seq_io import print_message
import sys

quality_score_encoding = {
    'phred33': {'!': 0, '\"': 1, '#': 2, '$': 3, '%': 4, '&': 5, '\'': 6,
        '(': 7, ')': 8, '*': 9, '+': 10, ',': 11, '-': 12, '.': 13, '/': 14,
        '0': 15, '1': 16, '2': 17, '3': 18, '4': 19, '5': 20, '6': 21, 
        '7': 22, '8': 23, '9': 24, ':': 25, ';': 26, '<': 27, '=': 28, 
        '>': 29, '?': 30, '@': 31, 'A': 32, 'B': 33, 'C': 34, 'D': 35, 
        'E': 36, 'F': 37, 'G': 38, 'H': 39, 'I': 40, 'J': 41, 'K': 42},
    'phred64': {'@': 0, 'A': 1, 'B': 2, 'C': 3, 'D': 4, 'E': 5, 'F': 6, 
        'G': 7, 'H': 8, 'I': 9, 'J': 10, 'K': 11, 'L': 12, 'M': 13, 'N': 14, 
        'O': 15, 'P': 16, 'Q': 17, 'R': 18, 'S': 19, 'T': 20, 'U': 21, 
        'V': 22, 'W': 23, 'X': 24, 'Y': 25, 'Z': 26, '[': 27, '\\': 28, 
        ']': 29, '^': 30, '_': 31, '`': 32, 'a': 33, 'b': 34, 'c': 35, 
        'd': 36, 'e': 37, 'f': 38, 'g': 39, 'h': 40, 'i': 41, 'j': 42}
    }

def translate_quality(base_qual, encoding='phred33'):
    """
    Translate ASCII characters to quality scores
    """
    try:
        base_score = quality_score_encoding[encoding][base_qual]
    except KeyError:
        print_message("error: wrong quality score encoding provided", 
            destination=sys.stderr)
    else:
        return base_score

def adaptive_trim(qual, trim_info, encoding='phred33'):
    """
    Uses sliding windows along with quality and length thresholds to determine \
    when the quality is sufficiently low to trim the 3'-end of reads and when \
    the quality is sufficiently high enough to trim the 5'-end of read
    (inspired by the trimmer sickle).
    """
    window_size = abs(trim_info[0])
    threshold = abs(trim_info[1])
    start = 0
    seqlen = end = len(qual)
    if seqlen == 0:
        return (start, end)

    if type(window_size) == type(int()):
        if window_size > seqlen or window_size == 1 or window_size == 0:
            step_size = seqlen
        else:
            step_size = window_size
    elif type(window_size) == type(float()):
        window_len = int(window_size * length)
        step_size = window_len if window_len > 1 else 2

    prev_scores = []
    found_start = False
    for position in range(start, end, step_size):
        frame = qual[position: position + step_size]
        framelen = len(frame)

        scores = [translate_quality(i, encoding) for i in frame]
        average = sum(scores) / framelen

        # find the start position by searching until the average > threshold
        if not found_start:
            if average > threshold:
                found_start = True
                # check to see if bases immediately before current frame are
                # above the threshold
                prev_scores.reverse()
                for score in prev_scores:
                    if score > threshold:
                        start -= 1
                    else:
                        break

                # average is lower than the threshold, but first few bases may
                # still be good quality
                if start == position:
                    for score in scores:
                        if score < threshold:
                            start += 1
                        else:
                            break
            else:
                start += framelen
        else:
            # now find the end position by searching until average < threshold
            if average < threshold:
                end = position
                # determine trim position by checking previous scores first
                prev_scores.reverse()
                for score in prev_scores:
                    if score < threshold:
                        end -= 1
                    else:
                        break

                # otherwise check scores of current frame and cut when it falls
                # below the threshold
                if end == position:
                    for score in scores:
                        if score >= threshold:
                            end += 1
                        else:
                            break

                return (start, end)

        prev_scores = scores

    # if no trimming required, return original start and end position
    return (start, end)

def trim_leading(qual, threshold, encoding='phred33'):
    """
    Trim low quality bases from the 5'-end of the sequence
    """
    threshold = abs(threshold)

    position = 0
    for position, base_qual in enumerate(qual):
        base_score = translate_quality(base_qual, encoding)
        if base_score >= threshold:
            break

    start = position
    return (start, len(qual))

def trim_trailing(qual, threshold, encoding='phred33'):
    """
    Trim low quality bases from the 3'-end of the sequence.
    """
    length = len(qual)
    threshold = abs(threshold)

    iterable = [-i for i in range(1, length + 1)]
    for position in iterable:
        base_score = translate_quality(qual[position], encoding)
        if base_score >= threshold:
            break
        else:
            length -= 1

    end = length
    return (0, end)

def trunc_n(seq):
    """
    Truncate sequence at first ambiguous base.
    """
    try:
        end = seq.index('N')
    except ValueError:
        end = len(seq)

    return (0, end)
