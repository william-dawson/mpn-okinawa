"""
MP1 integration test.

MP1 = one-body V operator contribution to correlation energy.
Result should be the HF exchange energy for a single V operator.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pdaggerq
from wick.helper import WickHelper


def pdaggerq_mp1():
    pq = pdaggerq.pq_helper("fermi")
    pq.set_left_operators([['1']])
    pq.add_operator_product(1.0, ['v'])
    pq.simplify()
    return pq.strings()


def ours_mp1():
    # Match pdaggerq reference path exactly: no extra linked-diagram filtering.
    wh = WickHelper(filter_unlinked=False, project_energy_subspace=False)
    wh.add_operator_product(1.0, ['v'])
    wh.simplify()
    return wh.strings()


def test_mp1():
    ref  = sorted(pdaggerq_mp1())
    ours = sorted(ours_mp1())
    assert ref == ours, f"\npdaggerq: {ref}\nours:     {ours}"

    print("PASS: MP1")


if __name__ == '__main__':
    test_mp1()
