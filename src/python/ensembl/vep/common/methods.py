
from collections import deque
import logging
import os

def trim_sequences(
    ref: str,
    alt: str,
    start: int = 0,
    end: int | None = None,
    empty_to_dash: bool = False,
    end_first: bool = False,
    strand: int = 1,
) -> tuple[str, str, int, int, bool]:
    """
    Trim common sequence from a ref/alt allele pair to give minimal alleles.

    Trims matching bases from the start then end (or end then start if
    end_first=True), adjusting coordinates accordingly.

    Args:
        ref: Reference allele sequence.
        alt: Alternate allele sequence.
        start: Start coordinate (default 0).
        end: End coordinate (default: start + len(ref) - 1).
        empty_to_dash: If True, replace empty alleles with "-".
        end_first: If True, trim from end before start.
        strand: 1 for forward, -1 for reverse strand.

    Returns:
        Tuple of (new_ref, new_alt, new_start, new_end, changed).

    Raises:
        ValueError: If ref or alt is None.
    """
    if ref is None or alt is None:
        raise ValueError("Missing reference or alternate sequence")

    if end is None:
        end = start + (len(ref) - 1)

    changed = False

    if end_first:
        # trim from right
        while ref and alt and ref[-1] == alt[-1]:
            ref = ref[:-1]
            alt = alt[:-1]
            if strand == -1:
                start += 1
            else:
                end -= 1
            changed = True

        # trim from left
        while ref and alt and ref[0] == alt[0]:
            ref = ref[1:]
            alt = alt[1:]
            if strand == -1:
                end -= 1
            else:
                start += 1
            changed = True

    else:
        # trim from left
        while ref and alt and ref[0] == alt[0]:
            ref = ref[1:]
            alt = alt[1:]
            if strand == -1:
                end -= 1
            else:
                start += 1
            changed = True

        # trim from right
        while ref and alt and ref[-1] == alt[-1]:
            ref = ref[:-1]
            alt = alt[:-1]
            if strand == -1:
                start += 1
            else:
                end -= 1
            changed = True

    if empty_to_dash:
        if ref == "":
            ref = "-"
        if alt == "":
            alt = "-"

    return (ref, alt, start, end, changed)

def shift_sequence(
    seq: str,
    post_seq: str,
    pre_seq: str,
    start: int,
    end: int,
    reverse: bool = False,
) -> tuple[int, str, int, int]:
    """
    Calculate the maximum possible 3' (or 5') shift for an allele sequence.

    Slides the allele along the flanking sequence as long as the
    leading/trailing base matches the next flanking base, updating
    coordinates and the HGVS allele string accordingly.

    Args:
        seq: The allele sequence to shift.
        post_seq: 3' flanking sequence (downstream).
        pre_seq: 5' flanking sequence (upstream).
        start: Variant start coordinate.
        end: Variant end coordinate.
        reverse: If True, shift in the 5' direction; otherwise 3'.

    Returns:
        Tuple of (shift_length, shifted_seq, new_start, new_end).
    """
    if seq is None or post_seq is None or pre_seq is None:
        raise ValueError("Missing sequence to shift")
    
    if seq == "":
        return (0, seq, start, end)

    seq_len = len(seq)
    shift_length = 0

    # Use deques for O(1) rotation instead of repeated string slicing
    seq_d = deque(seq)

    if reverse:
        pre_len = len(pre_seq)
        loop_limit = pre_len - seq_len + 1
        if loop_limit < 0:
            loop_limit = pre_len

        for n in range(1, loop_limit + 1):
            # Compare last char of allele with upstream flanking base
            if seq_d[-1] != pre_seq[pre_len - n]:
                break

            shift_length += 1
            start -= 1
            end -= 1

            # Rotate allele: move last char to front
            seq_d.rotate(1)
    else:
        post_len = len(post_seq)
        loop_limit = post_len - seq_len + 1
        if loop_limit < 0:
            loop_limit = post_len

        for n in range(0, loop_limit):
            # Compare first char of allele with downstream flanking base
            if seq_d[0] != post_seq[n]:
                break

            shift_length += 1
            start += 1
            end += 1

            # Rotate allele: move first char to end
            seq_d.rotate(-1)

    return (
        shift_length,
        "".join(seq_d),
        start,
        end,
    )

def setup_logging(level: str = "INFO"):
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(os.path.join(output_file, "_warnings.txt")),
        ]
    )