"""
Step 6c: Optional graph-based connected-diagram filter.

This pass builds a tensor connectivity graph for each term:
  - nodes: ERI tensors in the term
  - edge: two tensors share at least one index label

Only connected terms are kept.
"""

from __future__ import annotations
from typing import List
from .term import Term


def _is_connected(term: Term) -> bool:
    tensors = term.tensors
    n = len(tensors)
    if n <= 1:
        return True

    adj = {i: set() for i in range(n)}
    for i in range(n):
        a = set(tensors[i].indices)
        for j in range(i + 1, n):
            b = set(tensors[j].indices)
            if not a.isdisjoint(b):
                adj[i].add(j)
                adj[j].add(i)

    seen = {0}
    stack = [0]
    while stack:
        cur = stack.pop()
        for nxt in adj[cur]:
            if nxt not in seen:
                seen.add(nxt)
                stack.append(nxt)
    return len(seen) == n


def filter_connected_terms_graph(terms: List[Term]) -> List[Term]:
    """Keep only connected terms using tensor-index graph connectivity."""
    return [t for t in terms if _is_connected(t)]

