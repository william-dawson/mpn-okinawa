"""
Step 6a: Linked-Cluster Theorem Filter

Remove unlinked (disconnected) diagram contributions.

In MBPT, the linked-cluster theorem states that unlinked diagram contributions
to the energy cancel exactly. A diagram is LINKED if all interaction vertices
(tensors) are connected via shared indices.

Algorithm:
  - Build a connectivity graph where nodes = tensors, edges = shared indices
  - Check if all tensors form one connected component
  - Drop terms with multiple disconnected pieces

Example disconnected terms (REMOVED):
  - <i,j||i,j> * <k,l||k,l>          (two separate bubbles)
  - <i,j||a,b> * <k,l||k,l>          ((i,j,a,b) separate from (k,l))
  - <i,j||i,j> * <a,b||a,b> * <c,d||c,d>  (three separate pieces)

Example connected terms (KEPT):
  - <i,j||a,b> * <a,b||i,j>          (i,j,a,b all connected)
  - <i,j||a,b> * <a,k||b,l> * <k,l||i,j>  (all connected through indices)
"""

from __future__ import annotations
from typing import List
from .term import Term


def _tensors_share_index(t1, t2) -> bool:
    """Check if two tensors share at least one index."""
    return not set(t1.indices).isdisjoint(t2.indices)


def _is_connected(term: Term) -> bool:
    """Check if all tensors in a term are connected via shared indices."""
    tensors = term.tensors
    if len(tensors) <= 1:
        return True

    # BFS connectivity check: can we reach all tensors from tensor 0?
    visited = {0}
    stack = [0]

    while stack:
        i = stack.pop()
        for j in range(len(tensors)):
            if j in visited:
                continue
            if _tensors_share_index(tensors[i], tensors[j]):
                visited.add(j)
                stack.append(j)

    return len(visited) == len(tensors)


def filter_unlinked_diagrams(terms: List[Term]) -> List[Term]:
    """Keep only linked (connected) diagrams."""
    result = []
    for term in terms:
        if _is_connected(term):
            result.append(term)
    return result
