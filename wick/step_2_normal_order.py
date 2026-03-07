"""
Step 2: Normal Ordering via Wick's Theorem (Fermi Vacuum)

Algorithm
---------
Repeatedly scan each term for the leftmost pair of adjacent operators
(ops[i], ops[i+1]) that is out of Fermi-vacuum normal order, i.e. a
FV annihilator immediately followed by a FV creator.

When such a pair is found at positions i, i+1, two new terms are generated:

1. Swap term:
   Move ops[i+1] to position i and ops[i] to position i+1.
   Multiply sign by -1  (anticommutation: {A,B} = AB + BA -> AB = -BA + {A,B}).

2. Contraction term (only when the pair can contract):
   A FV annihilator '-X' and FV creator '+Y' can contract iff they are in
   the SAME orbital space (both occ or both virt), producing delta(X,Y).
   The operators at positions i and i+1 are removed; a Delta(X.label, Y.label)
   is appended.  Sign is unchanged.

   Contraction check (true vacuum perspective):
     FV annihilator of occ  = a†_i  (hole annihilator)
     FV creator of occ      = a_i   (hole creator)
     { a†_i, a_j } = delta_ij  -> contraction when BOTH OCC
     FV annihilator of virt = a_a   (particle annihilator)
     FV creator of virt     = a†_a  (particle creator)
     { a_a, a†_b } = delta_ab  -> contraction when BOTH VIRT

This is iterated (BFS over the tree of generated terms) until all terms
are in normal order.  Fully-contracted terms (operators list empty) and
non-zero-operator normal-ordered terms are both collected.
"""

from __future__ import annotations
from typing import List
from .term import Term, is_occ, is_virt, Delta


def _first_out_of_order(operators) -> int:
    """Return index i of leftmost FV annihilator before a FV creator, or -1."""
    for i in range(len(operators) - 1):
        if not operators[i].dagger_fermi and operators[i + 1].dagger_fermi:
            return i
    return -1


def _swap_one(term: Term) -> List[Term]:
    """
    Find the first out-of-order adjacent pair and return the 1 or 2 new terms.
    Returns the original term (unchanged) if already normal ordered.
    """
    ops = term.operators
    i   = _first_out_of_order(ops)

    if i == -1:
        return [term]  # already normal ordered

    ann = ops[i]      # FV annihilator
    cre = ops[i + 1]  # FV creator

    # 1. Swap term
    swap_ops  = ops[:i] + [cre, ann] + ops[i + 2:]
    swap_term = Term(
        term.factor,
        -term.sign,
        term.tensors,
        swap_ops,
        list(term.deltas),
    )

    results = [swap_term]

    # 2. Contraction term (same space only)
    both_occ  = is_occ(ann.label)  and is_occ(cre.label)
    both_virt = is_virt(ann.label) and is_virt(cre.label)

    if both_occ or both_virt:
        contr_ops  = ops[:i] + ops[i + 2:]
        contr_deltas = list(term.deltas) + [Delta(ann.label, cre.label)]
        contr_term = Term(
            term.factor,
            term.sign,
            term.tensors,
            contr_ops,
            contr_deltas,
        )
        results.append(contr_term)

    return results


def normal_order_fermi_vacuum(terms: List[Term]) -> List[Term]:
    """
    Bring all terms to Fermi-vacuum normal order using Wick's theorem.
    Returns ALL resulting terms (contracted and non-contracted alike).
    """
    ordered = []
    queue   = list(terms)

    while queue:
        term = queue.pop(0)
        if term.is_normal_ordered_fermi():
            ordered.append(term)
        else:
            queue.extend(_swap_one(term))

    return ordered
