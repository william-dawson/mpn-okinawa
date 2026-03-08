import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from wick.helper import WickHelper


@pytest.mark.parametrize(
    "n_ops,expected_parity,expected_default",
    [
        (1, 1, 1),
        (2, 2, 1),
        (3, 5, 3),
    ],
)
def test_mp_term_counts_match_current_cli_behavior(n_ops, expected_parity, expected_default):
    wh_parity = WickHelper(filter_unlinked=False)
    wh_parity.add_operator_product(1.0, ['v'] * n_ops)
    wh_parity.simplify()
    assert len(wh_parity.strings()) == expected_parity

    wh_default = WickHelper()
    wh_default.add_operator_product(1.0, ['v'] * n_ops)
    wh_default.simplify()
    assert len(wh_default.strings()) == expected_default


def test_denominators_align_with_output_terms():
    wh = WickHelper(filter_unlinked=False)
    wh.add_operator_product(1.0, ['v', 'v'])
    wh.simplify()

    rows = wh.strings()
    denoms = wh.denominators()
    assert len(rows) == len(denoms)
    assert len(rows) > 0
    assert all(isinstance(d, str) and len(d) > 0 for d in denoms)
