"""
Step 5: Label Canonicalization

Rename all summed (internal) index labels to conventional names:
  occupied  summed: i, j, k, l, m, n  (in order of first appearance)
  virtual   summed: a, b, c, d, e, f  (in order of first appearance)

"Summed" = not a free external label (for our vacuum-to-vacuum case, ALL
remaining labels are summed).

Labels are renamed in the order they appear scanning left to right across
tensors (then deltas, then operators) to give a deterministic result.
"""

from __future__ import annotations
from typing import List
from .term import Term, is_occ, is_virt


_OCC_NAMES  = list('ijklmn')
_VIRT_NAMES = list('abcdef')


def _collect_labels_in_order(term: Term) -> List[str]:
    """Return labels in the order they appear in the term (no duplicates)."""
    seen  = []
    order = []

    def visit(lbl: str):
        if lbl not in seen:
            seen.append(lbl)
            order.append(lbl)

    for t in term.tensors:
        for idx in t.indices:
            visit(idx)
    for d in term.deltas:
        visit(d.i)
        visit(d.j)
    for op in term.operators:
        visit(op.label)

    return order


def canonicalize_labels(terms: List[Term]) -> List[Term]:
    """Rename internal labels to i, j, k, ... and a, b, c, ... per term."""
    result = []
    for term in terms:
        labels = _collect_labels_in_order(term)

        occ_idx  = 0
        virt_idx = 0
        mapping  = {}

        for lbl in labels:
            if lbl in mapping:
                continue
            if is_occ(lbl):
                if occ_idx >= len(_OCC_NAMES):
                    raise ValueError(f"Too many occupied labels (max {len(_OCC_NAMES)})")
                mapping[lbl] = _OCC_NAMES[occ_idx]
                occ_idx += 1
            elif is_virt(lbl):
                if virt_idx >= len(_VIRT_NAMES):
                    raise ValueError(f"Too many virtual labels (max {len(_VIRT_NAMES)})")
                mapping[lbl] = _VIRT_NAMES[virt_idx]
                virt_idx += 1
            # else: leave unknown labels unchanged (shouldn't happen post-expand)

        result.append(term.renamed(mapping))

    return result


def canonicalize_labels_pdaggerq_style(terms: List[Term]) -> List[Term]:
    """
    Mirror pdaggerq use_conventional_labels behavior.

    Replace o<n> labels in increasing n order with first free occupied
    conventional labels, then replace v<n> labels similarly for virtuals.
    """
    occ_in = [f"o{i}" for i in range(30)]
    virt_in = [f"v{i}" for i in range(30)]
    occ_out = list("ijklmn")
    virt_out = list("abcdef")

    out = []
    for term in terms:
        used = set(_collect_labels_in_order(term))
        mapping = {}

        for src in occ_in:
            if src not in used:
                continue
            for dst in occ_out:
                if dst not in used and dst not in mapping.values():
                    mapping[src] = dst
                    used.add(dst)
                    break

        for src in virt_in:
            if src not in used:
                continue
            for dst in virt_out:
                if dst not in used and dst not in mapping.values():
                    mapping[src] = dst
                    used.add(dst)
                    break

        out.append(term.renamed(mapping))
    return out
