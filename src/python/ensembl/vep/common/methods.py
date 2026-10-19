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
import re
from typing import Any, Callable, Literal, TypeVar
from attrs import Attribute as AttrsAttribute

from ensembl.vep.common.enums import Strand, ClassSOTerm
from ensembl.vep.common.constants import CLASS_SO_TERM_ABBRV_TO_ENUM

logger = logging.getLogger(__name__)

T = TypeVar("T")


def is_resolved(instance: object, name: str) -> bool:
    """Return whether a slot-backed cached property is already resolved without triggering computation.

    Args:
        instance (object): Object that owns the cached property.
        name (str): Protected cached-property name (for example, "_cds").

    Returns:
        bool: True if the value is already present (seeded or previously computed),
            False if reading it would still trigger computation.

    Raises:
        AttributeError: If name is not a descriptor on the instance type hierarchy.
    """
    for cls in type(instance).__mro__:
        descriptor = cls.__dict__.get(name)
        if descriptor is not None and hasattr(descriptor, "__get__"):
            try:
                # Use descriptor.__get__ directly: normal access (including getattr)
                # would invoke __getattr__ and compute/cache the value on a miss.
                descriptor.__get__(instance, type(instance))  # pylint: disable=unnecessary-dunder-call
                return True
            except AttributeError:
                return False
    raise AttributeError(f"{name!r} is not a descriptor on {type(instance).__name__}")


def invalidate(instance: object, *invalidate_names: str) -> None:
    """Clear the named cached properties of a class object when the attributes are set.

    Args:
        instance (object): The object whose cached properties should be cleared.
        *invalidate_names (str): The names of the cached properties to invalidate.

    Returns:
        None
    """
    for name in invalidate_names:
        if not name.startswith("_"):
            logger.debug(
                f"Invalidation only works for protected attributes (starting with _), but got '{name}'"
            )
        with suppress(AttributeError):
            object.__delattr__(instance, f"{name}")


def invalidates(*cached_names: str) -> Callable:
    """Returns an attrs on_setattr hook that clears named cached properties.

    Args:
        *cached_names (str): The names of the cached properties to invalidate.

    Returns:
        Callable: An attrs on_setattr hook function.
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
            * `name` is the attribute name to reset
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
            name = attr_pair[0]
            default_value = attr_pair[1]
            if not name.startswith("_"):
                logger.debug(
                    "reset_attrs should only be used for protected attributes (starting with _),"
                    f" but got '{name}'."
                )
            else:
                object.__setattr__(instance, name, default_value)

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
    strand: Strand = Strand.FORWARD,
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
        raise ValueError("Missing reference or alternate sequence for trimming")

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
    strand: Strand = Strand.FORWARD,
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
    """

    changed = False

    if strand == Strand.UNKNOWN:
        logger.warning("Strand is unknown, cannot trim alleles with confidence. Returning untrimmed alleles.")
        return (ref, alt, start, end, changed)

    if trim_at == "left":
        while ref and alt and ref[0] == alt[0]:
            ref = ref[1:]
            alt = alt[1:]
            if strand == Strand.REVERSE:
                end -= 1
            else:
                start += 1
            changed = True
    else:
        while ref and alt and ref[-1] == alt[-1]:
            ref = ref[:-1]
            alt = alt[:-1]
            if strand == Strand.REVERSE:
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

    Raises:
        ValueError: If seq, post_seq, or pre_seq is None.
    """
    if seq is None or post_seq is None or pre_seq is None:
        raise ValueError("Missing sequence to shift")

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


def parse_region_string(region: str) -> tuple[str, int | None, int | None]:
    """Split a region string into its sequence-region name and optional bounds.

    The inverse of format_region_string: "1:100-200" -> ("1", 100, 200), "1:100" -> ("1", 100,
    None), "1" -> ("1", None, None). A name may itself contain ":" (as HLA contigs do), so only
    the last ":" separates the coordinates, and only when what follows parses as bounds.

    Args:
        region (str): A region string in the format "seq_region_name[:start[-end]]".

    Returns:
        tuple[str, int | None, int | None]: The name, start and end.

    Raises:
        ValueError: If region is empty, or its coordinates are not positive integers.
    """

    if not region:
        raise ValueError("region must not be empty")

    name, separator, coordinates = region.rpartition(":")
    if not separator:
        return region, None, None

    start_text, _, end_text = coordinates.partition("-")
    try:
        start = int(start_text)
        end = int(end_text) if end_text else None
    except ValueError as invalid:
        raise ValueError(f"region coordinates must be integers: {region!r}") from invalid

    if start < 1 or (end is not None and end < 1):
        raise ValueError("regions are 1-based and must be positive")

    return name, start, end


def format_region_string(
    region: str | None = None,
    *,
    seq_region_name: str | None = None,
    start: int | None = None,
    end: int | None = None,
) -> str:
    """Format a region string into "seq_region_name:start-end"
        from either an existing region or seq_region_name/start/end values.

    Args:
        region: A region string in the format "seq_region_name:start-end".
        seq_region_name: Sequence ID (e.g. chromosome name).
        start: Start coordinate (1-based).
        end: End coordinate (1-based).

    Returns:
        :A region string in the format "seq_region_name:start-end".
    """

    if region is not None:
        if any(value is not None for value in (seq_region_name, start, end)):
            raise ValueError("pass either region or seq_region_name/start/end, not both")

        return region

    if seq_region_name is None:
        raise ValueError("region requires either region or seq_region_name")

    if start is None:
        if end is not None:
            raise ValueError("end coordinate provided without start coordinate")
        return seq_region_name

    if end is None:
        return f"{seq_region_name}:{start}"

    if start < 1 or end < 1:
        raise ValueError("regions are 1-based and must be positive")

    return f"{seq_region_name}:{start}-{end}"


def infer_variant_so_term(allele_string: str, ref_correct: bool = False) -> ClassSOTerm | None:
    """
    Return the SO term for the class of the variation.

    Args:
        allele_string: The allele string for the variation.
        ref_correct: Whether the first allele is known to be the reference.

    Returns:
        The SO term string for the variation.
    """
    # default to sequence_alteration
    so_term = ClassSOTerm.SEQUENCE_ALTERATION

    # A/T, A/T/G
    if re.match(r"^[A-Z](/[A-Z])+$", allele_string):
        so_term = ClassSOTerm.SNV

    # (CAG)8/(CAG)9
    elif re.search(r"\)\d+", allele_string):
        so_term = ClassSOTerm.TANDEM_REPEAT

    else:
        alleles = allele_string.split("/")

        if len(alleles) > 1:
            ref = alleles[0]
            alts = alleles[1:]

            if ref == "-":
                if len(alts) == 1 and re.search(r"DEL", alts[0]):
                    # -/(LARGEDELETION)
                    so_term = ClassSOTerm.DELETION if ref_correct else ClassSOTerm.INDEL

                if not any(not re.match(r"^[A-Z]+$|INS", alt) for alt in alts):
                    # -/ATT, -/(LARGEINSERTION)
                    so_term = ClassSOTerm.INSERTION if ref_correct else ClassSOTerm.INDEL

            elif re.match(r"^[A-Z]+$", ref):
                if not any(not re.match(r"^-$|DEL", alt) for alt in alts):
                    # A/-, A/(LARGEDELETION)
                    so_term = ClassSOTerm.DELETION if ref_correct else ClassSOTerm.INDEL

                elif re.match(r"^[A-Z]+(/[A-Z]+)+$", allele_string):
                    # AA/TT => substitution, AA/TTTT => indel
                    same_size = all(len(alt) == len(ref) for alt in alts)
                    so_term = ClassSOTerm.SUBSTITUTION if same_size else ClassSOTerm.INDEL

            elif re.search(r"DEL", ref):
                if not any(alt != "-" for alt in alts):
                    # (LARGEDELETION)/-, (2345 BP DELETION)/-
                    so_term = ClassSOTerm.DELETION if ref_correct else ClassSOTerm.INDEL

        elif len(alleles) == 1:
            if re.search(r"INS", alleles[0]):
                # (LARGEINSERTION)
                so_term = ClassSOTerm.INSERTION if ref_correct else ClassSOTerm.INDEL
            elif re.search(r"DEL", alleles[0]):
                # (308 BP DELETION)
                so_term = ClassSOTerm.DELETION if ref_correct else ClassSOTerm.INDEL

    return so_term


def infer_symbolic_variant_so_term(allele_string: str, svtype: str | None = None) -> ClassSOTerm | None:
    """Infer the SO term for a symbolic structural variant.

    This mirrors the symbolic variant normalization performed by Perl VEP's
    ``get_SO_term``. A VCF-compliant symbolic allele takes precedence over
    ``svtype`` so that subtype information in values such as
    ``<INS:ME:ALU>`` is retained. ``svtype`` is used as a fallback for a
    missing or non-standard allele.

    Args:
        allele_string: Symbolic alternate allele or breakend string.
        svtype: Optional structural variant type from VCF INFO/SVTYPE.

    Returns:
        The inferred structural variant SO term, or ``None`` if the type is
        unsupported.
    """
    raw_type = allele_string.strip()
    if svtype and not re.match(r"^<?(?:DEL|INS|DUP|INV|CNV|CN=?[0-9]+)", raw_type, flags=re.IGNORECASE):
        raw_type = svtype.strip()

    # Normalize repeated variants, e.g. INS/INS or <INS>,<INS>.
    parts = re.split(r"[/,]", raw_type)
    normalized_parts = [part.strip().strip("<>").upper() for part in parts]
    if normalized_parts and all(part == normalized_parts[0] for part in normalized_parts):
        raw_type = parts[0].strip()

    mobile_element_match = re.search(r"(INS|DEL):(ME):?([A-Z0-9]+)?", raw_type, flags=re.IGNORECASE)
    if mobile_element_match:
        variant_class = mobile_element_match.group(1).upper()
        element = (mobile_element_match.group(3) or "ME").upper()
        element = "LINE1" if element == "L1" else element
        subtype = element if element in {"ALU", "HERV", "LINE1", "SVA"} else "ME"
        abbreviation = f"{variant_class}_{subtype}"
    elif re.search(r"DEL", raw_type, flags=re.IGNORECASE) and re.search(
        r"DUP", raw_type, flags=re.IGNORECASE
    ):
        abbreviation = "CNV"
    elif re.search(r"DUP:TANDEM", raw_type, flags=re.IGNORECASE):
        abbreviation = "TDUP"
    elif re.search(r"CNV:TR", raw_type, flags=re.IGNORECASE):
        abbreviation = "TREP"
    elif re.search(r"CN=?[0-9]", raw_type, flags=re.IGNORECASE):
        abbreviation = "CNV"
        if re.fullmatch(r"<?CN=?0>?", raw_type, flags=re.IGNORECASE):
            abbreviation = "DEL"
        elif re.fullmatch(r"<?CN=?2>?", raw_type, flags=re.IGNORECASE):
            abbreviation = "DUP"
    elif re.search(r"CNV", raw_type, flags=re.IGNORECASE):
        abbreviation = "CNV"
    elif re.search(r"[\[\]]|^\.|\.$", raw_type):
        abbreviation = "BND"
    elif raw_type.startswith("<") or raw_type.endswith(">"):
        abbreviation = re.sub(r"[<>]", "", raw_type).split(":", maxsplit=1)[0]
    else:
        abbreviation = raw_type

    return CLASS_SO_TERM_ABBRV_TO_ENUM.get(abbreviation.upper())


def reverse_complement(sequence: str) -> str:
    """Returns the reverse complement of a nucleotide sequence.

    Args:
        sequence (str): The sequence to reverse complement.

    Returns:
        str: The reverse complement of the sequence.
    """

    translation_table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(translation_table)[::-1]


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
