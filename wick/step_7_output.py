"""
Step 7: Output Formatting

Converts a list of fully-contracted Terms into the pdaggerq strings() format:
  [ [factor_str, tensor1_str, tensor2_str, ...], ... ]

Factor format matches pdaggerq: sign + decimal with minimum precision,
e.g. 0.25 -> '+0.2500', -0.125 -> '-0.1250'.
"""

from __future__ import annotations
from typing import List
from .term import Term


def _fmt_factor(value: float) -> str:
    """
    Format a float exactly as pdaggerq's minimum_precision() does.

    The algorithm (from pq_string.h):
      1. Work with the string of 10 * abs(value)
      2. Count decimal digits, tracking repeating digits
      3. Stop after 12 consecutive repeating digits
      4. Subtract the repeat count from precision if the last digit is '0'
      5. Enforce a minimum of 2 decimal places
    """
    sign    = '+' if value >= 0 else '-'
    abs_val = abs(value)
    scaled  = 10.0 * abs_val

    s = f'{scaled:.6f}'   # matches C++ std::fixed default precision

    precision               = 0
    decimal_point_found     = False
    is_repeated             = False
    last_digit              = ' '
    repeat_count            = 0

    for digit in s:
        is_repeated = (digit == last_digit)
        last_digit  = digit

        if digit == '.':
            decimal_point_found = True
        elif decimal_point_found and is_repeated:
            repeat_count += 1
            if repeat_count >= 12:
                break

        if not is_repeated:
            repeat_count = 0

        if decimal_point_found:
            precision += 1

    if precision >= repeat_count and last_digit == '0':
        precision -= repeat_count

    if precision < 2:
        precision = 2

    return sign + f'{abs_val:.{precision}f}'


def format_strings(terms: List[Term]) -> List[List[str]]:
    """Convert terms to pdaggerq-compatible list-of-list-of-string format."""
    result = []
    for term in terms:
        if not term.is_fully_contracted():
            continue  # skip uncontracted terms (should have been filtered already)
        row = [_fmt_factor(term.effective_factor)]
        row += [t.to_string() for t in term.tensors]
        result.append(row)
    return result
