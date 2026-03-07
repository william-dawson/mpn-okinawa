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
from .term import Term, Delta


def _count_in_tensors_and_ops(term: Term, label: str) -> int:
    count = 0
    for t in term.tensors:
        count += t.indices.count(label)
    for op in term.operators:
        if op.label == label:
            count += 1
    return count


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

    # Emulate pdaggerq gobble_deltas preference:
    # replace the first delta label if it appears in integrals/operators,
    # otherwise replace the second one.
    if _count_in_tensors_and_ops(term, i) > 0:
        keep, replace = j, i
    elif _count_in_tensors_and_ops(term, j) > 0:
        keep, replace = i, j
    else:
        keep, replace = j, i
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
