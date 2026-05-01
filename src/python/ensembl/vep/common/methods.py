# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Common functions."""

from collections import deque
from contextlib import suppress
import logging
import os
import sys
from typing import Any, Callable, Literal, TypeVar
from attrs import Attribute as AttrsAttribute

from .enums import KnownStrand

logger = logging.getLogger(__name__)

T = TypeVar("T")


def invalidate(instance: object, *invalidate_names: str) -> None:
    """Clears the named cached properties of a class object when the attributes are set.

    Args:
        instance: The class object to invalidate attributes of.
        *invalidate_names: The names of the cached properties to invalidate.
    """
    for name in invalidate_names:
        if not name.startswith("__"):
            logger.debug(
                f"Invalidation only works for private attributes (starting with __), but got '{name}'"
            )
        with suppress(AttributeError):
            object.__delattr__(instance, f"_{type(instance).__name__}{name}")


def invalidates(*cached_names: str) -> Callable:
    """
    Returns an attrs on_setattr hook that clears the named cached property whenever the attribute is set.

    Args:
        *cached_names: The names of the cached properties to invalidate.

    Returns:
        Callable: An on_setattr hook function.
    """

    def hook(instance: object, _modified_attr: AttrsAttribute, new_value: T) -> T:
        invalidate(instance, *cached_names)

        # attrs requires the hook to return the value
        return new_value

    return hook


def resets_attrs(*attrs: tuple[str, Any]) -> Callable:
    """
    Returns an attrs on_setattr hook that resets the named attributes to a default value.

    Args:
        *attrs: list of tuples (name, default_value) where
                 * `name` is the attribute name to resent
                 * `default_value` is the value to reset the attribute to on each call

    Returns:
        Callable: An on_setattr hook function.
    """

    def hook(instance: object, _modified_attr: AttrsAttribute, new_value: T) -> T:
        # nothing to do if the value hasn't changed
        if object.__getattribute__(instance, _modified_attr.name) == new_value:
            return new_value

        # Otherwise, reset the attributes defined in attrs
        for attr_pair in attrs:
            short_name = attr_pair[0]
            long_name = f"_{type(instance).__name__}{short_name}"
            default_value = attr_pair[1]
            if not short_name.startswith("__"):
                logger.debug(
                    "reset_attrs should only be used for private attributes (starting with __)"
                    ", but got '%s'.",
                    short_name,
                )
            else:
                object.__setattr__(instance, long_name, default_value)

        # attrs requires the hook to return the value on success
        return new_value

    return hook


def trim_sequences(
    ref: str,
    alt: str,
    *,
    start: int = 0,
    end: int | None = None,
    empty_to_dash: bool = False,
    end_first: bool = False,
    strand: KnownStrand = KnownStrand.FORWARD,
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
    """

    end = end or (start + (len(ref) - 1))

    changed = False

    if end_first:
        # trim from right
        ref, alt, start, end, right_changed = _trim_common_sequence(
            ref, alt, start=start, end=end, trim_at="right", strand=strand
        )
        changed = right_changed

        # trim from left
        ref, alt, start, end, left_changed = _trim_common_sequence(
            ref, alt, start=start, end=end, trim_at="left", strand=strand
        )
        changed = changed or left_changed

    else:
        # trim from left
        ref, alt, start, end, left_changed = _trim_common_sequence(
            ref, alt, start=start, end=end, trim_at="left", strand=strand
        )
        changed = left_changed

        # trim from right
        ref, alt, start, end, right_changed = _trim_common_sequence(
            ref, alt, start=start, end=end, trim_at="right", strand=strand
        )
        changed = changed or right_changed

    if empty_to_dash:
        ref = ref or "-"
        alt = alt or "-"

    return (ref, alt, start, end, changed)


def _trim_common_sequence(
    ref: str,
    alt: str,
    *,
    start: int = 0,
    end: int,
    trim_at: Literal["left", "right"],
    strand: KnownStrand = KnownStrand.FORWARD,
) -> tuple[str, str, int, int, bool]:
    """
    Trim common sequence from start or end of a ref/alt allele pair.

    Trims matching bases from the start or end, adjusting coordinates accordingly.

    Args:
        ref: Reference allele sequence.
        alt: Alternate allele sequence.
        start: Start coordinate (default 0).
        end: End coordinate (default: start + len(ref) - 1).
        trim_at: What side to trim from: 'left' or 'right'.
        strand: 1 for forward, -1 for reverse strand.

    Returns:
        Tuple of (new_ref, new_alt, new_start, new_end, changed).

    Raises:
        ValueError: If ref or alt is None.
    """

    changed = False

    if trim_at == "left":
        while ref and alt and ref[0] == alt[0]:
            ref = ref[1:]
            alt = alt[1:]
            if strand == KnownStrand.REVERSE:
                end -= 1
            else:
                start += 1
            changed = True
    else:
        while ref and alt and ref[-1] == alt[-1]:
            ref = ref[:-1]
            alt = alt[:-1]
            if strand == KnownStrand.REVERSE:
                start += 1
            else:
                end -= 1
            changed = True

    return (ref, alt, start, end, changed)


def shift_sequence(
    seq: str,
    *,
    post_seq: str,
    pre_seq: str,
    start: int,
    end: int,
    reverse: bool = False,
) -> tuple[int, str, int, int, bool]:
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
        Tuple of (shift_length, shifted_seq, new_start, new_end, changed).
    """

    if seq == "":
        return (0, seq, start, end, False)

    seq_len = len(seq)
    shift_length = 0

    # Use deques for O(1) rotation instead of repeated string slicing
    seq_d = deque(seq)

    changed = False

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

            changed = True
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

            changed = True

    return (shift_length, "".join(seq_d), start, end, changed)


def setup_logging(level: str = "INFO") -> None:
    """Set up logging configuration.

    Args:
        level (str): The logging level to use (default: "INFO").
    """
    output_file = "DUMMY_OUTPUT_FILE"
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(os.path.join(output_file, "_warnings.txt")),
        ],
    )
