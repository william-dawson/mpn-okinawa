"""
Parse pdaggerq strings() output into Term objects so we can feed them
to our pipeline steps.

pdaggerq format: [[factor_str, tensor1_str, tensor2_str, ...], ...]
  factor_str:  '+0.2500' or '-0.1250' etc.
  tensor_str:  '<p,q||r,s>'  for two-body
               'f(p,q)'      for one-body (future)
"""

from __future__ import annotations
import re
from typing import List
from .term import Term, Tensor

_ERI_RE = re.compile(r'^<([^,]+),([^|]+)\|\|([^,]+),([^>]+)>$')


def _parse_tensor_str(s: str) -> Tensor:
    m = _ERI_RE.match(s.strip())
    if m:
        return Tensor('g', [m.group(1), m.group(2), m.group(3), m.group(4)])
    raise ValueError(f"Cannot parse tensor string: {s!r}")


def pdaggerq_strings_to_terms(strings: List[List[str]]) -> List[Term]:
    """Convert pdaggerq strings() output to a list of Terms (no operators)."""
    terms = []
    for row in strings:
        factor_str = row[0]
        sign   = -1 if factor_str[0] == '-' else 1
        factor = abs(float(factor_str))
        tensors = [_parse_tensor_str(s) for s in row[1:]]
        terms.append(Term(factor, sign, tensors, operators=[]))
    return terms
