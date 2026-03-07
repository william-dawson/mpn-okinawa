"""
Step 1: General Label Expansion

pdaggerq expands operators with general indices (p, q, r, s, ...) into
occupied + virtual cases before applying Wick's theorem.  We do the same.

A term with n general labels expands into 2^n terms, one for each assignment
of each general label to either occupied (fresh occ label) or virtual
(fresh virt label).

Fresh label naming:
  occupied : o1, o2, o3, ...   (later canonicalized to i, j, k, ...)
  virtual  : v1, v2, v3, ...   (later canonicalized to a, b, c, ...)

The counters are per-call so that labels across different input terms
do not collide before canonicalization.
"""

from __future__ import annotations
from itertools import product
from typing import List

from .term import Term, is_general


_OCC_PREFIX  = 'o'
_VIRT_PREFIX = 'v'


def _general_labels_in_term(term: Term) -> List[str]:
    """Collect all unique general labels appearing in operators and tensors."""
    seen = []
    for op in term.operators:
        if is_general(op.label) and op.label not in seen:
            seen.append(op.label)
    for t in term.tensors:
        for idx in t.indices:
            if is_general(idx) and idx not in seen:
                seen.append(idx)
    return seen


def expand_general_labels(terms: List[Term], _counter: list = None) -> List[Term]:
    """
    Expand every general-label term into occ/virt-specific variants.

    A shared mutable counter list is used so that fresh labels are unique
    across all terms in a single simplify() call.
    """
    if _counter is None:
        _counter = [0]  # [next_fresh_index]

    result = []
    for term in terms:
        gen_labels = _general_labels_in_term(term)
        if not gen_labels:
            result.append(term)
            continue

        # Expand: for each general label assign 'occ' or 'virt'
        for assignment in product(('occ', 'virt'), repeat=len(gen_labels)):
            mapping = {}
            for label, space in zip(gen_labels, assignment):
                _counter[0] += 1
                n = _counter[0]
                if space == 'occ':
                    mapping[label] = f'{_OCC_PREFIX}{n}'
                else:
                    mapping[label] = f'{_VIRT_PREFIX}{n}'

            result.append(term.renamed(mapping))

    return result
