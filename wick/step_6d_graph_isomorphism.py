"""
Step 6d: Graph-isomorphism-based term consolidation.

This pass groups algebraically equivalent tensor products by a graph key
that is invariant to dummy-label naming and tensor ordering, then sums
their coefficients.
"""

from __future__ import annotations

from itertools import permutations
from typing import List

from .term import Term, Tensor, is_occ, is_virt

_EPSILON = 1e-12


def _canonicalize_tensor_internal(tensor: Tensor) -> tuple[Tensor, int]:
    """Canonicalize antisymmetric ERI slots within bra/ket pairs."""
    if tensor.name != "g" or len(tensor.indices) != 4:
        return tensor.copy(), 1

    p, q, r, s = tensor.indices
    sign = 1
    if p > q:
        p, q = q, p
        sign *= -1
    if r > s:
        r, s = s, r
        sign *= -1
    return Tensor("g", [p, q, r, s]), sign


def _canonicalize_term_internal(term: Term) -> tuple[Term, int]:
    tensors = []
    sign = 1
    for t in term.tensors:
        ct, s = _canonicalize_tensor_internal(t)
        tensors.append(ct)
        sign *= s
    return Term(term.factor, term.sign, tensors, []), sign


def _next_occ_name(n: int) -> str:
    base = "ijklmn"
    return base[n] if n < len(base) else f"i{n}"


def _next_virt_name(n: int) -> str:
    base = "abcdef"
    return base[n] if n < len(base) else f"a{n}"


def _canonical_graph_key(term: Term) -> tuple:
    """
    Exact canonical key under:
      - permutation of tensor factors
      - renaming of dummy occupied/virtual labels
    """
    n = len(term.tensors)
    if n == 0:
        return tuple()

    best = None
    for perm in permutations(range(n)):
        occ_map = {}
        virt_map = {}
        gen_map = {}
        occ_n = 0
        virt_n = 0
        gen_n = 0
        key_tensors = []

        for ti in perm:
            t = term.tensors[ti]
            mapped = []
            for lbl in t.indices:
                if is_occ(lbl):
                    if lbl not in occ_map:
                        occ_map[lbl] = _next_occ_name(occ_n)
                        occ_n += 1
                    mapped.append(occ_map[lbl])
                elif is_virt(lbl):
                    if lbl not in virt_map:
                        virt_map[lbl] = _next_virt_name(virt_n)
                        virt_n += 1
                    mapped.append(virt_map[lbl])
                else:
                    if lbl not in gen_map:
                        gen_map[lbl] = f"g{gen_n}"
                        gen_n += 1
                    mapped.append(gen_map[lbl])
            key_tensors.append((t.name, tuple(mapped)))

        key = tuple(key_tensors)
        if best is None or key < best:
            best = key

    return best


def _normalized_layout_key(term: Term) -> tuple[tuple[str, tuple[str, ...]], ...]:
    """
    Produce a deterministic representative layout by renaming labels in
    first-appearance order to conventional occ/virt names.
    """
    occ_names = iter("ijklmn")
    virt_names = iter("abcdef")
    mapping = {}

    def norm(lbl: str) -> str:
        if lbl in mapping:
            return mapping[lbl]
        if is_occ(lbl):
            mapping[lbl] = next(occ_names)
        elif is_virt(lbl):
            mapping[lbl] = next(virt_names)
        else:
            mapping[lbl] = lbl
        return mapping[lbl]

    rows = []
    for t in term.tensors:
        rows.append((t.name, tuple(norm(x) for x in t.indices)))
    rows.sort()
    return tuple(rows)


def cancel_terms_graph_isomorphic(terms: List[Term]) -> List[Term]:
    """
    Consolidate terms by graph isomorphism and sum coefficients.
    """
    buckets = {}

    for term in terms:
        canon_term, canon_sign = _canonicalize_term_internal(term)
        eff = term.effective_factor * canon_sign
        key = _canonical_graph_key(canon_term)
        rep_layout = _normalized_layout_key(canon_term)

        if key not in buckets:
            buckets[key] = {"sum": eff, "rep": rep_layout}
        else:
            buckets[key]["sum"] += eff
            if rep_layout < buckets[key]["rep"]:
                buckets[key]["rep"] = rep_layout

    out = []
    for data in buckets.values():
        total = data["sum"]
        if abs(total) <= _EPSILON:
            continue
        tensors = [Tensor(name, list(indices)) for name, indices in data["rep"]]
        out.append(Term(abs(total), 1 if total > 0 else -1, tensors, []))
    return out
