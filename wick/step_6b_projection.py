"""
Step 6b: Energy-subspace projection.

For vacuum-to-vacuum energy expressions built from two-body interactions,
terms with an odd occupied/virtual balance in the final contracted tensors
correspond to single-excitation-like couplings and should be removed.

This mirrors the operator-portion elimination behavior used by pdaggerq's
simplify path for these MP-style energy expressions.
"""

from __future__ import annotations
from typing import List
from .term import Term, is_occ, is_virt


def _label_balance(term: Term) -> tuple[int, int]:
    labels = set()
    for tensor in term.tensors:
        labels.update(tensor.indices)
    n_occ = sum(1 for lbl in labels if is_occ(lbl))
    n_virt = sum(1 for lbl in labels if is_virt(lbl))
    return n_occ, n_virt


def filter_energy_subspace(terms: List[Term]) -> List[Term]:
    """Drop terms with odd occupied/virtual label balance."""
    out = []
    for term in terms:
        n_occ, n_virt = _label_balance(term)
        if (n_occ % 2) != 0 or (n_virt % 2) != 0:
            continue
        out.append(term)
    return out

