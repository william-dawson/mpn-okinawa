"""
Step 7 isolated test: cancel terms.

cancel_terms() sums factors of terms with identical tensor structures and
drops any term whose total factor is (effectively) zero.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from wick.term import Term, Tensor
from wick.step_6_cleanup import cancel_terms


def test_canonicalize_and_merge():
    # After canonicalization, <i,j||a,b> and <a,b||i,j> become the same: <a,b||i,j>
    # So these two 0.25 factors should merge to 0.50
    t1 = Term(0.25, 1, [Tensor('g', ['i','j','a','b'])], [])
    t2 = Term(0.25, 1, [Tensor('g', ['a','b','i','j'])], [])
    result = cancel_terms([t1, t2])
    assert len(result) == 1
    assert abs(result[0].effective_factor - 0.50) < 1e-12
    print("PASS: canonicalization merges <i,j||a,b> and <a,b||i,j>")


def test_identical_terms_summed():
    # Two identical terms: 0.25 + 0.25 = 0.50
    t1 = Term(0.25, 1, [Tensor('g', ['i','j','a','b'])], [])
    t2 = Term(0.25, 1, [Tensor('g', ['i','j','a','b'])], [])
    result = cancel_terms([t1, t2])
    assert len(result) == 1
    assert abs(result[0].effective_factor - 0.50) < 1e-12
    print("PASS: identical terms summed to 0.50")


def test_opposite_terms_cancel():
    t1 = Term(0.25,  1, [Tensor('g', ['i','j','a','b'])], [])
    t2 = Term(0.25, -1, [Tensor('g', ['i','j','a','b'])], [])
    result = cancel_terms([t1, t2])
    assert result == []
    print("PASS: opposite-sign terms cancel to zero")


def test_partial_cancellation():
    # 0.50 - 0.25 = 0.25
    t1 = Term(0.50,  1, [Tensor('g', ['i','j','a','b'])], [])
    t2 = Term(0.25, -1, [Tensor('g', ['i','j','a','b'])], [])
    result = cancel_terms([t1, t2])
    assert len(result) == 1
    assert abs(result[0].effective_factor - 0.25) < 1e-12
    print("PASS: partial cancellation (0.50 - 0.25 = 0.25)")


def test_sign_of_result_negative():
    # 0.25 - 0.50 = -0.25
    t1 = Term(0.25,  1, [Tensor('g', ['i','j','a','b'])], [])
    t2 = Term(0.50, -1, [Tensor('g', ['i','j','a','b'])], [])
    result = cancel_terms([t1, t2])
    assert len(result) == 1
    assert result[0].sign == -1
    assert abs(result[0].factor - 0.25) < 1e-12
    print("PASS: negative result has sign=-1")


def test_multi_tensor_canonical_merge():
    # After canonicalization and sorting, both tensor orders become the same
    t1 = Term(0.25, 1, [Tensor('g', ['a','b','i','j']), Tensor('g', ['i','j','a','b'])], [])
    t2 = Term(0.25, 1, [Tensor('g', ['a','b','i','j']), Tensor('g', ['i','j','a','b'])], [])
    t3 = Term(0.25, 1, [Tensor('g', ['i','j','a','b']), Tensor('g', ['a','b','i','j'])], [])
    result = cancel_terms([t1, t2, t3])
    # All three canonicalize to [<a,b||i,j>, <a,b||i,j>], so merge to one term
    assert len(result) == 1
    assert abs(result[0].effective_factor - 0.75) < 1e-12
    print("PASS: multi-tensor canonical merge (0.25*3=0.75)")


def test_empty_list():
    assert cancel_terms([]) == []
    print("PASS: empty list")


if __name__ == '__main__':
    test_canonicalize_and_merge()
    test_identical_terms_summed()
    test_opposite_terms_cancel()
    test_partial_cancellation()
    test_sign_of_result_negative()
    test_multi_tensor_canonical_merge()
    test_empty_list()
