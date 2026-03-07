"""
Step 3: Contraction Filter

For the Fermi vacuum with vacuum bra <0| and ket |0>, only fully-contracted
strings (no remaining operators) survive.  All other terms vanish.
"""

from __future__ import annotations
from typing import List
from .term import Term


def filter_fully_contracted(terms: List[Term]) -> List[Term]:
    """Keep only terms with no remaining uncontracted operators."""
    return [t for t in terms if t.is_fully_contracted()]
