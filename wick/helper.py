"""
WickHelper: main entry point mirroring the pq_helper interface.

Input API
---------
    wh = WickHelper()
    wh.add_term(factor, tensors=['g(p,q,r,s)', ...], ops=['+p', '+q', '-s', '-r', ...])
    wh.simplify()
    result = wh.strings()   # list of [factor_str, tensor_str, ...]

Operator format
---------------
  '+label' -> a†_label  (true vacuum creator)
  '-label' -> a_label   (true vacuum annihilator)

Index conventions
-----------------
  i, j, k, l, m, n  -> occupied
  a, b, c, d, e, f  -> virtual
  anything else      -> general (expanded to occ + virt before Wick's theorem)

Pipeline (in order)
-------------------
  1. expand_general_labels   - expand p,q,r,... into occ+virt cases
  2. normal_order_fermi_vacuum - Wick's theorem
  3. filter_fully_contracted  - drop uncontracted terms
  4. apply_deltas             - substitute and remove delta functions
  5. canonicalize_labels      - rename to i,j,k.../a,b,c...
  6. cancel_terms             - sum identical terms, drop zeros
  7. format_strings           - convert to pdaggerq string format
"""

from __future__ import annotations
import re
from typing import List

from .term                  import Term, Tensor, Operator, is_occ, is_virt
from .step_1_expand         import expand_general_labels
from .step_2_normal_order   import normal_order_fermi_vacuum
from .step_3_filter         import filter_fully_contracted
from .step_4_deltas         import apply_deltas
from .step_5_labels         import canonicalize_labels
from .step_6_cleanup        import cancel_terms
from .step_6a_linked_cluster import filter_unlinked_diagrams
from .step_6b_projection    import filter_energy_subspace
from .step_7_output         import format_strings
from .step_8_denominator    import extract_denominator


# ── input parsing helpers ────────────────────────────────────────────────────

_TENSOR_RE = re.compile(r'^(\w+)\(([^)]+)\)$')


def _parse_tensor(s: str) -> Tensor:
    """Parse 'g(p,q,r,s)' into Tensor('g', ['p','q','r','s'])."""
    s = s.strip()
    m = _TENSOR_RE.match(s)
    if not m:
        raise ValueError(f"Cannot parse tensor string: {s!r}")
    name    = m.group(1)
    indices = [x.strip() for x in m.group(2).split(',')]
    return Tensor(name, indices)


def _parse_op(s: str) -> Operator:
    """Parse '+p' or '-p' into Operator('p', dagger=True/False)."""
    s = s.strip()
    if s[0] == '+':
        return Operator(s[1:], dagger=True)
    if s[0] == '-':
        return Operator(s[1:], dagger=False)
    raise ValueError(f"Operator string must start with '+' or '-', got: {s!r}")


# ── main class ───────────────────────────────────────────────────────────────

class WickHelper:

    def __init__(self, filter_unlinked: bool = True, project_energy_subspace: bool = True):
        self._input_terms: List[Term] = []
        self._result:      List[Term] = []
        self._filter_unlinked = filter_unlinked
        self._project_energy_subspace = project_energy_subspace

    def add_term(
        self,
        factor:  float,
        tensors: List[str],
        ops:     List[str],
    ) -> None:
        """
        Add one second-quantized input term.

        Parameters
        ----------
        factor  : numerical coefficient
        tensors : list of tensor strings, e.g. ['g(p,q,r,s)', 'g(t,u,v,w)']
        ops     : list of operator strings, e.g. ['+p', '+q', '-s', '-r']
                  The ops list is the PRODUCT of ALL operators for this term
                  (including contributions from multiple tensors if needed).
        """
        parsed_tensors = [_parse_tensor(t) for t in tensors]
        parsed_ops     = [_parse_op(o)     for o in ops]
        self._input_terms.append(Term(factor, 1, parsed_tensors, parsed_ops))

    def simplify(self) -> None:
        """Run the full pipeline."""
        terms = list(self._input_terms)

        # 1. Expand general labels (p,q,r,... -> occ/virt pairs)
        counter = [0]
        terms = expand_general_labels(terms, counter)

        # 2. Normal ordering (Wick's theorem, Fermi vacuum)
        terms = normal_order_fermi_vacuum(terms)

        # 3. Keep only fully contracted terms
        terms = filter_fully_contracted(terms)

        # 4. Apply delta functions
        terms = apply_deltas(terms)

        # 5. Canonicalize index labels
        terms = canonicalize_labels(terms)

        # 6. Cancel duplicate terms / sum factors
        terms = cancel_terms(terms)

        # 6b. Project out odd occupied/virtual balance terms (MP energy subspace)
        if self._project_energy_subspace:
            terms = filter_energy_subspace(terms)

        # 6a. Optional: filter unlinked (disconnected) diagrams via linked-cluster theorem
        if self._filter_unlinked:
            terms = filter_unlinked_diagrams(terms)

        self._result = terms

    def strings(self) -> List[List[str]]:
        """Return results in pdaggerq strings() format."""
        return format_strings(self._result)

    def denominators(self) -> List[str]:
        """Return energy denominators for each term in canonical form."""
        return [extract_denominator(term) for term in self._result]

    def terms_with_denominators(self) -> List[tuple]:
        """Return (strings_row, denominator) pairs for each result term."""
        result = []
        strings = self.strings()
        denoms = self.denominators()
        for string_row, denom in zip(strings, denoms):
            result.append((string_row, denom))
        return result

    def clear(self) -> None:
        self._input_terms = []
        self._result      = []
