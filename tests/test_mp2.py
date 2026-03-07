"""
MP2 integration test.

ours_mp2() starts as a direct call to pdaggerq (trivially passes).
As each pipeline step is verified and swapped in, this function is updated.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pdaggerq
from wick.helper import WickHelper


def pdaggerq_mp2():
    pq = pdaggerq.pq_helper("fermi")
    pq.set_left_operators([['1']])
    pq.add_operator_product(1.0, ['v', 'v'])
    pq.simplify()
    return pq.strings()


def ours_mp2():
    # Match pdaggerq reference path exactly: no extra linked-diagram filtering.
    wh = WickHelper(filter_unlinked=False)
    wh.add_term(
        1/16,
        ['g(p,q,r,s)', 'g(t,u,v,w)'],
        ['+p','+q','-r','-s', '+t','+u','-v','-w'],
    )
    wh.simplify()
    return wh.strings()


def test_mp2():
    ref  = sorted(pdaggerq_mp2())
    ours = sorted(ours_mp2())
    assert ref == ours, f"\npdaggerq: {ref}\nours:     {ours}"
    print("PASS: MP2")


if __name__ == '__main__':
    test_mp2()
