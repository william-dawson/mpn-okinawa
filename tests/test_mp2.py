"""
MP2 integration test.

ours_mp2() starts as a direct call to pdaggerq (trivially passes).
As each pipeline step is verified and swapped in, this function is updated.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pdaggerq
from wick.helper import WickHelper
from wick.term import is_occ, is_virt
from wick.step_6_cleanup import _canonicalize_integral_indices
import re

_ERI_RE = re.compile(r'^<([^,]+),([^|]+)\|\|([^,]+),([^>]+)>$')


def _normalize_rows(rows):
    bucket = {}

    for row in rows:
        factor = float(row[0])
        tensors = []
        sign = 1
        for s in row[1:]:
            canon, sgn = _canonicalize_integral_indices(s)
            sign *= sgn
            m = _ERI_RE.match(canon)
            tensors.append([m.group(1), m.group(2), m.group(3), m.group(4)])

        tensors.sort()
        mapping = {}
        occ_names = iter("ijklmn")
        virt_names = iter("abcdef")
        key_tensors = []
        for t in tensors:
            kt = []
            for lbl in t:
                if lbl not in mapping:
                    if is_occ(lbl):
                        mapping[lbl] = next(occ_names)
                    elif is_virt(lbl):
                        mapping[lbl] = next(virt_names)
                    else:
                        mapping[lbl] = lbl
                kt.append(mapping[lbl])
            key_tensors.append(tuple(kt))
        key = tuple(key_tensors)
        bucket[key] = bucket.get(key, 0.0) + factor * sign

    return {k: round(v, 12) for k, v in bucket.items() if abs(v) > 1e-12}


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
    ref = _normalize_rows(pdaggerq_mp2())
    ours = _normalize_rows(ours_mp2())
    assert ref == ours, f"\npdaggerq(norm): {ref}\nours(norm):     {ours}"
    print("PASS: MP2")


if __name__ == '__main__':
    test_mp2()
