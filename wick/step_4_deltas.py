"""
Step 4: Delta Elimination

Each Delta(i, j) in a term means delta_{ij}.  We substitute one label for
the other everywhere (tensors, remaining deltas) and remove the delta.

Substitution rule: replace the label that is a fresh internal label
(o1, v2, ...) with the other one, preferring to keep conventional labels
(i, j, k, ..., a, b, c, ...) in the output.

If both labels are fresh, keep the one with the smaller numeric suffix.
If both are conventional, keep the one that comes first alphabetically
(pdaggerq keeps one and drops the other).

A delta delta(i, i) (same label on both sides) just gets removed with no
substitution.
"""

from __future__ import annotations
from typing import List
from .term import Term, Delta, is_occ, is_virt


def _is_fresh(label: str) -> bool:
    """Fresh labels have the form o<n> or v<n>."""
    return (label.startswith('o') or label.startswith('v')) and label[1:].isdigit()


def _prefer_keep(a: str, b: str) -> tuple:
    """Return (keep, replace) choosing which label to eliminate."""
    a_fresh = _is_fresh(a)
    b_fresh = _is_fresh(b)

    if a_fresh and not b_fresh:
        return b, a   # keep conventional b, replace fresh a
    if b_fresh and not a_fresh:
        return a, b   # keep conventional a, replace fresh b
    if a_fresh and b_fresh:
        # both fresh: keep lower numeric index
        na = int(a[1:])
        nb = int(b[1:])
        if na <= nb:
            return a, b
        return b, a
    # both conventional: keep alphabetically first
    if a <= b:
        return a, b
    return b, a


def _apply_one_delta(term: Term) -> Term:
    """Apply the first delta in term.deltas and return the updated term."""
    if not term.deltas:
        return term

    delta = term.deltas[0]
    remaining_deltas = list(term.deltas[1:])

    i, j = delta.i, delta.j

    if i == j:
        # delta(x, x) = 1, just remove it
        return Term(term.factor, term.sign, term.tensors, term.operators, remaining_deltas)

    keep, replace = _prefer_keep(i, j)
    mapping = {replace: keep}

    new_tensors = [t.renamed(mapping) for t in term.tensors]
    new_deltas  = [d.renamed(mapping) for d in remaining_deltas]
    new_ops     = [op.renamed(mapping) for op in term.operators]

    return Term(term.factor, term.sign, new_tensors, new_ops, new_deltas)


def apply_deltas(terms: List[Term]) -> List[Term]:
    """Iteratively apply all delta functions in each term."""
    result = []
    for term in terms:
        while term.deltas:
            term = _apply_one_delta(term)
        result.append(term)
    return result
