"""
Step 4 isolated test: label canonicalization.

Given terms with fresh expanded labels (o1, o2, v1, v2, ...), verify
canonicalize_labels() renames them to i,j,k,... and a,b,c,... in
first-appearance order across the tensors.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from wick.term import Term, Tensor
from wick.step_5_labels import canonicalize_labels


def test_occ_only():
    # g(o1,o2,o1,o2) -> g(i,j,i,j)
    t = Term(1.0, 1, [Tensor('g', ['o1','o2','o1','o2'])], [])
    result = canonicalize_labels([t])
    assert result[0].tensors[0].indices == ['i','j','i','j'], result[0].tensors[0].indices
    print("PASS: occ-only labels -> i,j")


def test_virt_occ():
    # g(v1,v2,o1,o2) -> g(a,b,i,j)
    t = Term(1.0, 1, [Tensor('g', ['v1','v2','o1','o2'])], [])
    result = canonicalize_labels([t])
    assert result[0].tensors[0].indices == ['a','b','i','j'], result[0].tensors[0].indices
    print("PASS: virt/occ labels -> a,b,i,j")


def test_two_tensors_shared_labels():
    # g(o1,o2,v1,v2) g(v1,v2,o1,o2) -> g(i,j,a,b) g(a,b,i,j)
    t = Term(1.0, 1,
             [Tensor('g', ['o1','o2','v1','v2']),
              Tensor('g', ['v1','v2','o1','o2'])], [])
    result = canonicalize_labels([t])
    assert result[0].tensors[0].indices == ['i','j','a','b'], result[0].tensors[0].indices
    assert result[0].tensors[1].indices == ['a','b','i','j'], result[0].tensors[1].indices
    print("PASS: two tensors with shared labels")


def test_first_appearance_order():
    # Virtual appears before occupied in the first tensor: g(v3,o5,v3,o5)
    # First label seen is v3 -> a, then o5 -> i
    t = Term(1.0, 1, [Tensor('g', ['v3','o5','v3','o5'])], [])
    result = canonicalize_labels([t])
    assert result[0].tensors[0].indices == ['a','i','a','i'], result[0].tensors[0].indices
    print("PASS: first-appearance order (virt before occ)")


def test_factor_sign_preserved():
    t = Term(0.25, -1, [Tensor('g', ['o3','o7','v2','v5'])], [])
    result = canonicalize_labels([t])
    assert result[0].factor == 0.25
    assert result[0].sign   == -1
    assert result[0].tensors[0].indices == ['i','j','a','b']
    print("PASS: factor and sign preserved")


def test_disconnected_term():
    # Two separate tensors with independent labels:
    # g(o1,o2,o1,o2) g(o3,o4,o3,o4) -> g(i,j,i,j) g(k,l,k,l)
    t = Term(0.25, 1,
             [Tensor('g', ['o1','o2','o1','o2']),
              Tensor('g', ['o3','o4','o3','o4'])], [])
    result = canonicalize_labels([t])
    assert result[0].tensors[0].indices == ['i','j','i','j'], result[0].tensors[0].indices
    assert result[0].tensors[1].indices == ['k','l','k','l'], result[0].tensors[1].indices
    print("PASS: disconnected term (independent labels)")


if __name__ == '__main__':
    test_occ_only()
    test_virt_occ()
    test_two_tensors_shared_labels()
    test_first_appearance_order()
    test_factor_sign_preserved()
    test_disconnected_term()
