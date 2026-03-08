"""
Step 4b: Reclassify fluctuation-potential one-body pieces.

pdaggerq introduces an intermediate occupied-repulsion object from the j1
piece of the fluctuation potential, then reclassifies it into an ERI by
introducing an extra occupied summation index:

    <?i||?j>  ->  <i,k||j,k>

This step mirrors that behavior for terms carrying tensors with
name == 'occ_repulsion'.
"""

from __future__ import annotations
from typing import List
from .term import Term, Tensor


_OCC_PREF = ["i", "j", "k", "l", "m", "n"]


def _all_labels(term: Term) -> set[str]:
    labels = set()
    for t in term.tensors:
        labels.update(t.indices)
    for d in term.deltas:
        labels.add(d.i)
        labels.add(d.j)
    for op in term.operators:
        labels.add(op.label)
    return labels


def _fresh_occ_label(used: set[str]) -> str:
    for lbl in _OCC_PREF:
        if lbl not in used:
            return lbl
    n = 0
    while True:
        lbl = f"o{n}"
        if lbl not in used:
            return lbl
        n += 1


def reclassify_occ_repulsion(terms: List[Term]) -> List[Term]:
    """Convert occ_repulsion tensors into equivalent ERI tensors."""
    out = []
    for term in terms:
        used = _all_labels(term)
        new_tensors = []
        for t in term.tensors:
            if t.name != "occ_repulsion":
                new_tensors.append(t.copy())
                continue
            i, j = t.indices
            k = _fresh_occ_label(used)
            used.add(k)
            new_tensors.append(Tensor("g", [i, k, j, k]))
        out.append(Term(term.factor, term.sign, new_tensors, [op.copy() for op in term.operators], [d.copy() for d in term.deltas]))
    return out

