from __future__ import division

from seq_qc.seq_io import print_error
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

def adaptive_trim(scores, trim_info):
    """
    Uses sliding window along with quality and length thresholds to determine \
    when the quality is sufficiently low to trim the 3'-end of reads and when \
    the quality is sufficiently high enough to trim the 5'-end of read. 
    Inspired by the trimmer sickle.
    """
    # Number of bases to remove from either end
    start, end = 0, 0
    seqlen = len(scores)
    if seqlen == 0:
        return (start, end)

    window_size = abs(trim_info[0])
    threshold = abs(trim_info[1])
    if type(window_size) == type(int()):
        if window_size > seqlen or window_size == 0:
            def_size = 0.1 * seqlen
            step_size = round(def_size) if def_size >= 1 else seqlen
        else:
            step_size = window_size
    elif type(window_size) == type(float()) and (0 < window_size < 1):
        window_len = round(window_size * seqlen)
        step_size = window_len if window_len > 1 else 2
    else:
        def_size = 0.1 * seqlen
        step_size = round(def_size) if def_size >= 1 else seqlen

    prev_scores = []
    found_start = False
    for position in range(start, seqlen, step_size):
        frame = scores[position: position + step_size]
        framelen = len(frame)
        average = sum(frame) / framelen

        # find the start position by searching until the average > threshold
        if not found_start:
            if average > threshold:
                found_start = True
                # check to see if bases immediately before current frame are
                # above the threshold too
                for score in prev_scores[::-1]:
                    if score > threshold:
                        start -= 1
                    else:
                        break

                # average is higher than the threshold, but the rise in quality
                # may have occurred towards the end of the frame
                if start == position:
                    for score in frame:
                        if score < threshold:
                            start += 1
                        else:
                            break
            else:
                start += framelen
        else:
            # now find the end position by searching until average < threshold
            if average < threshold:
                end = seqlen - position
                # determine trim position by checking previous scores first
                for score in prev_scores[::-1]:
                    if score < threshold:
                        end += 1
                    else:
                        break

                # otherwise check scores of current frame and cut when it falls
                # below the threshold
                if end == seqlen - position:
                    for score in frame:
                        if score >= threshold:
                            end -= 1
                        else:
                            break

                return (start, end)

        prev_scores = frame

    # if no trimming required, return original start and end position
    return (start, end)

def trim_leading(scores, threshold):
    """
    Trim low quality bases from the 5'-end of the sequence
    """
    threshold = abs(threshold)

    position = 0
    for position, basescore in enumerate(scores):
        if basescore >= threshold:
            break

    return (position, 0)

def trim_trailing(scores, threshold):
    """
    Trim low quality bases from the 3'-end of the sequence.
    """
    threshold = abs(threshold)

    position = 0
    for position, basescore in enumerate(scores[::-1]):
        if basescore >= threshold:
            break

    return (0, position)


def truncate_by_n(record):
    """
    Truncate reads to the first position containing an ambiguous base
    """
    try:
        nstart = record.sequence.index('N')
        record.sequence, record.quality = record.sequence[0: nstart], \
            record.quality[0: nstart]
        return record
    except ValueError:
        return record
