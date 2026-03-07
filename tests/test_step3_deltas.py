"""
Step 3 isolated test: delta elimination.

apply_deltas() iteratively substitutes delta_{ij} by replacing one label with
the other everywhere (tensors, remaining deltas, operators) and removes the delta.

Substitution rule:
  - fresh label (o<n>/v<n>) replaced by conventional (i,j,... / a,b,...) if present
  - both fresh: lower numeric index survives
  - both conventional: alphabetically first survives
  - delta(x, x): removed with no substitution
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from wick.term import Term, Tensor, Delta, Operator
from wick.step_4_deltas import apply_deltas


def test_trivial_delta():
    # delta(o1, o1) = 1 — just removed
    t = Term(0.25, 1, [Tensor('g', ['o1','o2','o1','o2'])], [], [Delta('o1','o1')])
    result = apply_deltas([t])
    assert result[0].deltas == []
    assert result[0].tensors[0].indices == ['o1','o2','o1','o2']
    print("PASS: trivial delta(x,x)")


def test_fresh_fresh_lower_survives():
    # delta(o1, o2) with g(o2, v1, o3, v2) -> replace o2 with o1
    t = Term(1.0, 1, [Tensor('g', ['o2','v1','o3','v2'])], [], [Delta('o1','o2')])
    result = apply_deltas([t])
    assert result[0].deltas == []
    assert result[0].tensors[0].indices == ['o1','v1','o3','v2']
    print("PASS: fresh+fresh: lower index (o1) survives")


def test_fresh_replaced_by_conventional():
    # delta(o1, i) -> replace o1 with i everywhere
    t = Term(1.0, 1, [Tensor('g', ['o1','v1','o2','v2'])], [], [Delta('o1','i')])
    result = apply_deltas([t])
    assert result[0].deltas == []
    assert result[0].tensors[0].indices == ['i','v1','o2','v2']
    print("PASS: fresh replaced by conventional")


def test_conventional_replaced_by_conventional():
    # delta(j, i) -> keep alphabetically first (i), replace j
    t = Term(1.0, 1, [Tensor('g', ['j','b','i','a'])], [], [Delta('j','i')])
    result = apply_deltas([t])
    assert result[0].deltas == []
    assert result[0].tensors[0].indices == ['i','b','i','a']
    print("PASS: conventional+conventional: alphabetically first survives")


def test_multiple_deltas_sequential():
    # delta(o1,o3) delta(o2,o4) with g(o1,o2,v1,v2) g(v1,v2,o3,o4)
    # After d(o1,o3): g(o1,o2,v1,v2) g(v1,v2,o1,o4), d(o2,o4)
    # After d(o2,o4): g(o1,o2,v1,v2) g(v1,v2,o1,o2)
    t = Term(0.25, 1,
             [Tensor('g', ['o1','o2','v1','v2']),
              Tensor('g', ['v1','v2','o3','o4'])],
             [],
             [Delta('o1','o3'), Delta('o2','o4')])
    result = apply_deltas([t])
    assert result[0].deltas == []
    assert result[0].tensors[0].indices == ['o1','o2','v1','v2']
    assert result[0].tensors[1].indices == ['v1','v2','o1','o2']
    print("PASS: two deltas applied sequentially")


def test_delta_propagates_through_remaining_deltas():
    # delta(o1,o2) delta(o2,o3) ->
    # first: replace o2 with o1 in remaining deltas and tensors
    #        -> delta(o1,o3), g with o2->o1
    # second: replace o3 with o1
    t = Term(1.0, 1,
             [Tensor('g', ['o2','v1','o3','v2'])],
             [],
             [Delta('o1','o2'), Delta('o2','o3')])
    result = apply_deltas([t])
    assert result[0].deltas == []
    # o2 -> o1 in tensor: g(o1, v1, o3, v2)
    # o2 -> o1 in second delta: Delta(o1, o3)
    # then apply Delta(o1,o3): o3 -> o1 in tensor: g(o1, v1, o1, v2)
    assert result[0].tensors[0].indices == ['o1','v1','o1','v2']
    print("PASS: delta propagates through remaining deltas")


def test_factor_sign_unchanged():
    t = Term(0.125, -1, [Tensor('g', ['o1','v1','o2','v2'])], [], [Delta('o1','o2')])
    result = apply_deltas([t])
    assert result[0].factor == 0.125
    assert result[0].sign == -1
    print("PASS: factor and sign unchanged by delta elimination")


def test_delta_in_operators():
    # delta(o1, o2) should also rename operators
    op = Operator('o1', True)
    t = Term(1.0, 1, [Tensor('g', ['o1','v1','o2','v2'])], [op], [Delta('o1','o2')])
    result = apply_deltas([t])
    assert result[0].deltas == []
    # o2 replaced by o1 everywhere (o1 < o2)
    assert result[0].operators[0].label == 'o1'
    print("PASS: delta renames operator labels")


if __name__ == '__main__':
    test_trivial_delta()
    test_fresh_fresh_lower_survives()
    test_fresh_replaced_by_conventional()
    test_conventional_replaced_by_conventional()
    test_multiple_deltas_sequential()
    test_delta_propagates_through_remaining_deltas()
    test_factor_sign_unchanged()
    test_delta_in_operators()
