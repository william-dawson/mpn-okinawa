"""
Step 6: Cleanup / Simplification with ERI Canonicalization

Two passes:
  1. Canonicalize each integral within each tensor (sort indices, apply sign).
  2. Collect terms by canonical form and sum their numerical factors.
  3. Drop any term whose total factor is (effectively) zero.

Canonicalization rules for <p,q||r,s>:
  - Sort bra (p,q) lexicographically, applying -1 if swapped
  - Sort ket (r,s) lexicographically, applying -1 if swapped
  - Do not swap bra/ket blocks
"""

from __future__ import annotations
from itertools import permutations
import re
from typing import List
from .term import Term, Tensor, is_occ, is_virt


_EPSILON = 1e-12


def _canonicalize_integral_indices(integral_str: str) -> tuple[str, int]:
    """
    Transform '<p,q||r,s>' into canonical form with sign.
    Returns (canonical_string, sign_change).
    """
    indices = re.findall(r'[a-z0-9]+', integral_str)
    if len(indices) != 4:
        return (integral_str, 1)

    p, q, r, s = indices
    sign = 1

    # Rule 1: sort bra (p,q)
    if p > q:
        p, q = q, p
        sign *= -1

    # Rule 2: sort ket (r,s)
    if r > s:
        r, s = s, r
        sign *= -1

    return (f"<{p},{q}||{r},{s}>", sign)


def _canonicalize_tensor(tensor: Tensor) -> tuple[Tensor, int]:
    """
    Canonicalize all 4-index integrals in a tensor.
    Returns (new_tensor, total_sign).
    """
    name = tensor.name
    indices = tensor.indices
    sign = 1

    if len(indices) == 4:
        integral_str = f"<{indices[0]},{indices[1]}||{indices[2]},{indices[3]}>"
        canon_str, s = _canonicalize_integral_indices(integral_str)
        sign *= s

        canon_indices = re.findall(r'[a-z0-9]+', canon_str)
        return (Tensor(name, canon_indices), sign)

    return (tensor, sign)


def _normalize_labels(
    tensors: List[tuple[str, tuple[str, ...]]]
) -> tuple[tuple[tuple[str, tuple[str, ...]], ...], int]:
    """
    Canonicalize dummy index names so equivalent summed-index terms merge.
    """
    occ_names = iter("ijklmn")
    virt_names = iter("abcdef")
    mapping = {}

    def map_label(lbl: str) -> str:
        if lbl in mapping:
            return mapping[lbl]
        if is_occ(lbl):
            mapping[lbl] = next(occ_names)
        elif is_virt(lbl):
            mapping[lbl] = next(virt_names)
        else:
            mapping[lbl] = lbl
        return mapping[lbl]

    normalized = []
    sign = 1
    for name, idxs in tensors:
        mapped = [map_label(x) for x in idxs]
        if len(mapped) == 4:
            t, s = _canonicalize_tensor(Tensor(name, mapped))
            sign *= s
            normalized.append((name, tuple(t.indices)))
        else:
            normalized.append((name, tuple(mapped)))
    normalized.sort()
    return tuple(normalized), sign


def _term_key(term: Term) -> tuple[tuple[tuple[str, tuple[str, ...]], ...], int]:
    """
    Canonical key with sign, including permutation of occupied/virtual dummy labels.
    """
    occ_labels = sorted({idx for t in term.tensors for idx in t.indices if is_occ(idx)})
    virt_labels = sorted({idx for t in term.tensors for idx in t.indices if is_virt(idx)})

    best_key = None
    best_sign = 1

    for occ_perm in permutations(occ_labels) if occ_labels else [()]:
        occ_map = dict(zip(occ_labels, occ_perm))
        for virt_perm in permutations(virt_labels) if virt_labels else [()]:
            rename_map = dict(occ_map)
            rename_map.update(dict(zip(virt_labels, virt_perm)))

            canon_tensors = []
            total_sign = 1
            for tensor in term.tensors:
                renamed = tensor.renamed(rename_map)
                canon_tensor, s = _canonicalize_tensor(renamed)
                canon_tensors.append((canon_tensor.name, tuple(canon_tensor.indices)))
                total_sign *= s

            canon_tensors.sort()
            key, s2 = _normalize_labels(canon_tensors)
            total_sign *= s2

            if best_key is None or key < best_key or (key == best_key and total_sign > best_sign):
                best_key = key
                best_sign = total_sign

    return best_key, best_sign


def cancel_terms(terms: List[Term]) -> List[Term]:
    """Canonicalize integrals, sum identical terms, drop zeros."""
    buckets = {}   # key -> running_factor

    for term in terms:
        key, canon_sign = _term_key(term)
        effective = term.effective_factor * canon_sign

        if key in buckets:
            buckets[key] += effective
        else:
            buckets[key] = effective

    result = []
    for key, total_factor in buckets.items():
        if abs(total_factor) > _EPSILON:
            tensors = [Tensor(name, list(indices)) for name, indices in key]
            new_term = Term(abs(total_factor), 1 if total_factor > 0 else -1, tensors, [])
            result.append(new_term)

    return result
