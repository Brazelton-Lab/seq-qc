from __future__ import division
from seq_io import print_error
import sys

def translate_quality(quals, encoding=33):
    """
    Translate ASCII characters to quality scores
    """
    valid_range = range(0, 43)
    qscores = [ord(i) - encoding for i in quals]
    for qscore in qscores:
        if qscore not in valid_range:
            print_error("error: wrong quality score encoding provided")
    return qscores

def adaptive_trim(quals, trim_info, encoding=33):
    """
    Uses sliding windows along with quality and length thresholds to determine \
    when the quality is sufficiently low to trim the 3'-end of reads and when \
    the quality is sufficiently high enough to trim the 5'-end of read
    (inspired by the trimmer sickle).
    """
    window_size = abs(trim_info[0])
    threshold = abs(trim_info[1])
    start = 0
    seqlen = end = len(quals)
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
        frame = quals[position: position + step_size]
        framelen = len(frame)

        scores = translate_quality(frame, encoding)
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
                # still be good qualsity
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

def trim_leading(quals, threshold, encoding=33):
    """
    Trim low quality bases from the 5'-end of the sequence
    """
    threshold = abs(threshold)

    position = 0
    for position, basescore in enumerate(translate_quality(quals, encoding)):
        if basescore >= threshold:
            break

    start = position
    return (start, len(quals))

def trim_trailing(quals, threshold, encoding=33):
    """
    Trim low quality bases from the 3'-end of the sequence.
    """
    end = len(quals)
    threshold = abs(threshold)

    scores = translate_quality(quals, encoding)
    scores.reverse()
    for basescore in scores:
        if basescore >= threshold:
            break
        else:
            end -= 1

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
