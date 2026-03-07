"""
Core data structures for second-quantized operator algebra.

Operator input convention
--------------------------
Operators are specified in TRUE VACUUM notation:
  '+p'  -> a†_p  (creation operator)
  '-p'  -> a_p   (annihilation operator)

The Fermi vacuum type (dagger_fermi) is derived automatically from the label:
  a†_virt  ->  FV creator     (dagger_fermi = True)
  a†_occ   ->  FV annihilator (dagger_fermi = False)  -- fills a hole
  a_virt   ->  FV annihilator (dagger_fermi = False)
  a_occ    ->  FV creator     (dagger_fermi = True)   -- creates a hole

Index label conventions
-----------------------
  Occupied : labels starting with i, j, k, l, m, n
  Virtual  : labels starting with a, b, c, d, e, f
  General  : anything else (p, q, r, s, t, u, v, w, ...)
             General labels must be expanded before dagger_fermi is accessed.
"""

from __future__ import annotations
from typing import List, Optional

# ── index classification ────────────────────────────────────────────────────

_OCC_CHARS  = frozenset('ijklmn')
_VIRT_CHARS = frozenset('abcdef')


def is_occ(label: str) -> bool:
    # conventional: i,j,k,l,m,n...  OR  fresh expanded: o1, o2, ...
    return label[0] in _OCC_CHARS or (label[0] == 'o' and label[1:].isdigit())

def is_virt(label: str) -> bool:
    # conventional: a,b,c,d,e,f...  OR  fresh expanded: v1, v2, ...
    return label[0] in _VIRT_CHARS or (label[0] == 'v' and label[1:].isdigit())

def is_general(label: str) -> bool:
    return not is_occ(label) and not is_virt(label)


# ── Operator ────────────────────────────────────────────────────────────────

class Operator:
    """A single fermionic creation/annihilation operator (true vacuum basis)."""

    __slots__ = ('label', 'dagger')

    def __init__(self, label: str, dagger: bool):
        self.label  = label   # index label, e.g. 'p', 'i1', 'a'
        self.dagger = dagger  # True = a†_label, False = a_label (TRUE vacuum)

    @property
    def dagger_fermi(self) -> bool:
        """Creator relative to Fermi vacuum?

        dagger_fermi = dagger XOR is_occ(label)
        """
        if is_occ(self.label):
            return not self.dagger
        if is_virt(self.label):
            return self.dagger
        raise ValueError(
            f"dagger_fermi undefined for general label '{self.label}'. "
            "Expand general labels first."
        )

    def renamed(self, mapping: dict) -> 'Operator':
        return Operator(mapping.get(self.label, self.label), self.dagger)

    def copy(self) -> 'Operator':
        return Operator(self.label, self.dagger)

    def __repr__(self) -> str:
        return f"{'+'if self.dagger else '-'}{self.label}"

    def __eq__(self, other) -> bool:
        return isinstance(other, Operator) and self.label == other.label and self.dagger == other.dagger


# ── Tensor ──────────────────────────────────────────────────────────────────

class Tensor:
    """A named tensor with index labels."""

    __slots__ = ('name', 'indices')

    def __init__(self, name: str, indices: List[str]):
        self.name    = name
        self.indices = list(indices)

    def renamed(self, mapping: dict) -> 'Tensor':
        return Tensor(self.name, [mapping.get(i, i) for i in self.indices])

    def copy(self) -> 'Tensor':
        return Tensor(self.name, list(self.indices))

    def to_string(self) -> str:
        """Format as pdaggerq-style string."""
        if self.name == 'g':
            p, q, r, s = self.indices
            return f'<{p},{q}||{r},{s}>'
        if self.name == 'f':
            p, q = self.indices
            return f'f({p},{q})'
        return f'{self.name}({",".join(self.indices)})'

    def __repr__(self) -> str:
        return self.to_string()

    def __eq__(self, other) -> bool:
        return isinstance(other, Tensor) and self.name == other.name and self.indices == other.indices


# ── Delta ───────────────────────────────────────────────────────────────────

class Delta:
    """Kronecker delta  δ_{ij}."""

    __slots__ = ('i', 'j')

    def __init__(self, i: str, j: str):
        self.i = i
        self.j = j

    def renamed(self, mapping: dict) -> 'Delta':
        return Delta(mapping.get(self.i, self.i), mapping.get(self.j, self.j))

    def copy(self) -> 'Delta':
        return Delta(self.i, self.j)

    def __repr__(self) -> str:
        return f'd({self.i},{self.j})'

    def __eq__(self, other) -> bool:
        return isinstance(other, Delta) and self.i == other.i and self.j == other.j


# ── Term ────────────────────────────────────────────────────────────────────

class Term:
    """A single term: sign * factor * Π tensors * Π deltas * Π operators."""

    def __init__(
        self,
        factor:    float,
        sign:      int,
        tensors:   List[Tensor],
        operators: List[Operator],
        deltas:    Optional[List[Delta]] = None,
    ):
        self.factor    = factor
        self.sign      = sign          # +1 or -1
        self.tensors   = list(tensors)
        self.operators = list(operators)
        self.deltas    = list(deltas) if deltas else []

    @property
    def effective_factor(self) -> float:
        return self.sign * self.factor

    # ── bookkeeping ─────────────────────────────────────────────────────────

    def copy(self) -> 'Term':
        return Term(
            self.factor,
            self.sign,
            [t.copy() for t in self.tensors],
            [op.copy() for op in self.operators],
            [d.copy() for d in self.deltas],
        )

    def renamed(self, mapping: dict) -> 'Term':
        return Term(
            self.factor,
            self.sign,
            [t.renamed(mapping) for t in self.tensors],
            [op.renamed(mapping) for op in self.operators],
            [d.renamed(mapping) for d in self.deltas],
        )

    # ── normal-order queries ─────────────────────────────────────────────────

    def is_normal_ordered_fermi(self) -> bool:
        """True iff all FV creators precede all FV annihilators."""
        seen_annihilator = False
        for op in self.operators:
            df = op.dagger_fermi
            if not df:
                seen_annihilator = True
            elif seen_annihilator:
                return False
        return True

    def is_fully_contracted(self) -> bool:
        return len(self.operators) == 0

    def has_general_labels(self) -> bool:
        for op in self.operators:
            if is_general(op.label):
                return True
        for t in self.tensors:
            for idx in t.indices:
                if is_general(idx):
                    return True
        return False

    def __repr__(self) -> str:
        parts = [f'{self.effective_factor:+.6f}']
        parts += [t.to_string() for t in self.tensors]
        parts += [repr(d) for d in self.deltas]
        parts += [repr(op) for op in self.operators]
        return ' * '.join(parts)
