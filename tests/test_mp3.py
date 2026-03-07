"""
MP3 integration test.

ours_mp3() starts as a direct call to pdaggerq (trivially passes).
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pdaggerq
from wick.helper import WickHelper


def pdaggerq_mp3():
    pq = pdaggerq.pq_helper("fermi")
    pq.set_left_operators([['1']])
    pq.add_operator_product(1.0, ['v', 'v', 'v'])
    pq.simplify()
    return pq.strings()


def ours_mp3():
    # Match pdaggerq reference path exactly: no extra linked-diagram filtering.
    wh = WickHelper(filter_unlinked=False)
    wh.add_term(
        1/64,
        ['g(p,q,r,s)', 'g(t,u,v,w)', 'g(x,y,z,o)'],
        ['+p','+q','-r','-s', '+t','+u','-v','-w', '+x','+y','-z','-o'],
    )
    wh.simplify()
    return wh.strings()


def test_mp3():
    ref  = sorted(pdaggerq_mp3())
    ours = sorted(ours_mp3())
    assert ref == ours, f"\npdaggerq: {ref}\nours:     {ours}"
    print("PASS: MP3")


if __name__ == '__main__':
    test_mp3()
