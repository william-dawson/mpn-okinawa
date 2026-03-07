"""
Step 8: Extract energy denominators from final terms.

After the full pipeline, each term represents a contribution to the correlation energy.
The denominator determines the orbital energy dependence: (sum of occupied energies) - (sum of virtual energies).

For MP2 with indices i,j (occupied) and a,b (virtual):
  denominator = (e_i + e_j) - (e_a + e_b) = e_i + e_j - e_a - e_b

This step extracts which indices participate (and with what sign).
"""

from __future__ import annotations
from typing import List, Set, Tuple
from .term import Term, is_occ, is_virt


def extract_denominator(term: Term) -> str:
    """
    Extract the energy denominator for a single term.

    Returns a string like "e_i + e_j - e_a - e_b" representing the denominator.
    Only unique indices are included (first appearance in left-to-right scan).
    """
    # Collect all indices appearing in tensors, in order of first appearance
    seen: Set[str] = set()
    occ_indices: List[str] = []
    virt_indices: List[str] = []

    for tensor in term.tensors:
        for idx in tensor.indices:
            if idx not in seen:
                seen.add(idx)
                if is_occ(idx):
                    occ_indices.append(idx)
                elif is_virt(idx):
                    virt_indices.append(idx)

    # Build denominator string: e_i + e_j - e_a - e_b
    parts = []

    for idx in occ_indices:
        parts.append(f"e_{idx}")

    for idx in virt_indices:
        parts.append(f"-e_{idx}")

    if not parts:
        return "1"  # No indices, denominator is 1

    # Join with proper formatting
    result = parts[0]
    for part in parts[1:]:
        if part.startswith('-'):
            result += f" {part}"
        else:
            result += f" + {part}"

    return result


def compute_denominators(terms: List[Term]) -> List[Tuple[Term, str]]:
    """
    Compute energy denominators for a list of terms.

    Returns list of (term, denominator_string) tuples.
    """
    return [(term, extract_denominator(term)) for term in terms]
