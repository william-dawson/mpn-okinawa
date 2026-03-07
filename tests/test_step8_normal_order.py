"""
Step 8 isolated test: Fermi-vacuum normal ordering (Wick's theorem).

Core rules being tested:
  - FV annihilator before FV creator is out of order
  - Out-of-order pair (A, B) generates two terms:
      swap:        sign flipped, operators swapped
      contraction: same sign, pair removed, Delta(ann.label, cre.label) added
        -- but ONLY when both operators are in the same orbital space
           (both occ OR both virt), i.e. when true-vacuum dagger flags differ
  - Cross-space out-of-order pair generates only the swap term (no contraction)

FV type recap (dagger_fermi = dagger XOR is_occ):
  +i (a†_i, occ) -> FV annihilator   -i (a_i, occ)  -> FV creator
  +a (a†_a, vir) -> FV creator        -a (a_a, vir)  -> FV annihilator
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from wick.term import Term, Tensor, Operator, Delta
from wick.step_2_normal_order import normal_order_fermi_vacuum, _swap_one
from wick.step_3_filter import filter_fully_contracted
from wick.step_4_deltas import apply_deltas


def _make_term(*op_strs, factor=1.0, sign=1, tensors=None):
    """Build a Term from operator strings like '+i', '-a', '+o1'."""
    ops = []
    for s in op_strs:
        dagger = (s[0] == '+')
        ops.append(Operator(s[1:], dagger))
    return Term(factor, sign, tensors or [], ops)


# ── _swap_one unit tests ────────────────────────────────────────────────────

def test_already_normal_ordered_no_ops():
    t = Term(1.0, 1, [Tensor('g', ['i','j','a','b'])], [])
    result = _swap_one(t)
    assert len(result) == 1
    assert result[0] is t
    print("PASS: no operators -> already normal ordered")


def test_already_normal_ordered_virt_creator_then_annihilator():
    # +a (FV cre) then -b (FV ann) -> normal ordered
    t = _make_term('+a', '-b')
    result = _swap_one(t)
    assert len(result) == 1
    assert result[0] is t
    print("PASS: +a -b -> already normal ordered")


def test_already_normal_ordered_occ_creator_then_annihilator():
    # -i (FV cre) then +j (FV ann) -> normal ordered
    t = _make_term('-i', '+j')
    result = _swap_one(t)
    assert len(result) == 1
    assert result[0] is t
    print("PASS: -i +j -> already normal ordered")


def test_occ_out_of_order_generates_swap_and_contraction():
    # +i (FV ann) then -j (FV cre): same space (occ), daggers differ -> swap + contraction
    t = _make_term('+i', '-j')
    result = _swap_one(t)
    assert len(result) == 2

    swap_t  = result[0]
    contr_t = result[1]

    # swap: sign flipped, operators reversed
    assert swap_t.sign == -1
    assert swap_t.operators[0].label == 'j' and not swap_t.operators[0].dagger  # -j first
    assert swap_t.operators[1].label == 'i' and swap_t.operators[1].dagger      # +i second
    assert swap_t.deltas == []

    # contraction: same sign, no operators, delta(i, j)
    assert contr_t.sign == 1
    assert contr_t.operators == []
    assert len(contr_t.deltas) == 1
    assert contr_t.deltas[0].i == 'i' and contr_t.deltas[0].j == 'j'
    print("PASS: +i -j -> swap (sign=-1) + contraction delta(i,j)")


def test_virt_out_of_order_generates_swap_and_contraction():
    # -a (FV ann) then +b (FV cre): same space (virt), daggers differ -> swap + contraction
    t = _make_term('-a', '+b')
    result = _swap_one(t)
    assert len(result) == 2

    swap_t  = result[0]
    contr_t = result[1]

    assert swap_t.sign == -1
    assert swap_t.operators[0].label == 'b' and swap_t.operators[0].dagger   # +b first
    assert swap_t.operators[1].label == 'a' and not swap_t.operators[1].dagger  # -a second

    assert contr_t.sign == 1
    assert contr_t.operators == []
    assert contr_t.deltas[0].i == 'a' and contr_t.deltas[0].j == 'b'
    print("PASS: -a +b -> swap (sign=-1) + contraction delta(a,b)")


def test_cross_space_generates_only_swap():
    # +i (FV ann, occ) then +a (FV cre, virt): daggers same (both True) -> swap only
    t = _make_term('+i', '+a')
    result = _swap_one(t)
    assert len(result) == 1
    swap_t = result[0]
    assert swap_t.sign == -1
    assert swap_t.operators[0].label == 'a'   # +a first
    assert swap_t.operators[1].label == 'i'   # +i second
    assert swap_t.deltas == []
    print("PASS: +i +a (cross-space) -> swap only, no contraction")


def test_cross_space_no_contraction_other_direction():
    # -a (FV ann, virt) then -i (FV cre, occ): both non-dagger -> swap only
    t = _make_term('-a', '-i')
    result = _swap_one(t)
    assert len(result) == 1
    assert result[0].sign == -1
    assert result[0].deltas == []
    print("PASS: -a -i (cross-space) -> swap only, no contraction")


def test_sign_propagates_from_parent():
    # parent has sign=-1; contraction should keep -1, swap should flip to +1
    t = _make_term('+i', '-j', sign=-1)
    result = _swap_one(t)
    swap_t, contr_t = result
    assert swap_t.sign == 1     # -1 * -1 = +1
    assert contr_t.sign == -1   # stays -1
    print("PASS: sign from parent propagates correctly")


def test_preceding_operators_preserved():
    # [+a, -b, +i, -j]: first two are normal ordered (+a FV cre, -b FV ann)
    # out-of-order pair is at positions 2,3: +i (FV ann) then -j (FV cre)
    t = _make_term('+a', '-b', '+i', '-j')
    result = _swap_one(t)
    assert len(result) == 2
    swap_t, contr_t = result
    # swap: [+a, -b, -j, +i]
    assert [o.label for o in swap_t.operators] == ['a', 'b', 'j', 'i']
    # contraction: [+a, -b] + delta(i,j)
    assert [o.label for o in contr_t.operators] == ['a', 'b']
    assert contr_t.deltas[0].i == 'i'
    print("PASS: operators before out-of-order pair are preserved")


# ── normal_order_fermi_vacuum integration-level tests ──────────────────────

def test_fully_contracted_vev_occ():
    # <0| a†_i a_i |0>_FV = 1
    # Input: [+i, -i] -> after NO + filter + delta -> effective_factor = +1
    t = _make_term('+i', '-i')
    no_terms  = normal_order_fermi_vacuum([t])
    contracted = filter_fully_contracted(no_terms)
    assert len(contracted) == 1
    final = apply_deltas(contracted)
    assert abs(final[0].effective_factor - 1.0) < 1e-12
    print("PASS: <0|a†_i a_i|0>_FV = 1")


def test_vev_virt_gives_zero():
    # <0| a†_a a_a |0>_FV = 0  (no virtual particles in FV)
    # Input: [+a, -a] is already normal ordered (+a FV cre, -a FV ann)
    # filter removes it since operators remain
    t = _make_term('+a', '-a')
    no_terms   = normal_order_fermi_vacuum([t])
    contracted = filter_fully_contracted(no_terms)
    assert contracted == []
    print("PASS: <0|a†_a a_a|0>_FV = 0")


def test_vev_occ_different_labels():
    # <0| a†_i a_j |0>_FV = delta_ij
    # After NO + filter: 1 term with delta(i,j)
    t = _make_term('+i', '-j')
    no_terms   = normal_order_fermi_vacuum([t])
    contracted = filter_fully_contracted(no_terms)
    assert len(contracted) == 1
    assert len(contracted[0].deltas) == 1
    d = contracted[0].deltas[0]
    assert {d.i, d.j} == {'i', 'j'}
    assert contracted[0].effective_factor > 0
    print("PASS: <0|a†_i a_j|0>_FV gives delta(i,j) term")


def test_vev_virt_different_labels():
    # <0| a_a a†_b |0>_FV = delta_ab
    t = _make_term('-a', '+b')
    no_terms   = normal_order_fermi_vacuum([t])
    contracted = filter_fully_contracted(no_terms)
    assert len(contracted) == 1
    assert len(contracted[0].deltas) == 1
    d = contracted[0].deltas[0]
    assert {d.i, d.j} == {'a', 'b'}
    assert contracted[0].effective_factor > 0
    print("PASS: <0|a_a a†_b|0>_FV gives delta(a,b) term")


def test_no_contraction_across_spaces():
    # <0| a†_i a†_a |0>_FV = 0  (can't contract occ with virt)
    t = _make_term('+i', '+a')
    no_terms   = normal_order_fermi_vacuum([t])
    contracted = filter_fully_contracted(no_terms)
    assert contracted == []
    print("PASS: <0|a†_i a†_a|0>_FV = 0 (no cross-space contraction)")


def test_factor_preserved_through_normal_order():
    t = _make_term('+i', '-j', factor=0.25, sign=-1)
    no_terms = normal_order_fermi_vacuum([t])
    contracted = filter_fully_contracted(no_terms)
    assert len(contracted) == 1
    # effective_factor should be -0.25 (sign kept through contraction)
    assert abs(contracted[0].effective_factor - (-0.25)) < 1e-12
    print("PASS: factor and sign preserved through normal ordering")


if __name__ == '__main__':
    test_already_normal_ordered_no_ops()
    test_already_normal_ordered_virt_creator_then_annihilator()
    test_already_normal_ordered_occ_creator_then_annihilator()
    test_occ_out_of_order_generates_swap_and_contraction()
    test_virt_out_of_order_generates_swap_and_contraction()
    test_cross_space_generates_only_swap()
    test_cross_space_no_contraction_other_direction()
    test_sign_propagates_from_parent()
    test_preceding_operators_preserved()
    test_fully_contracted_vev_occ()
    test_vev_virt_gives_zero()
    test_vev_occ_different_labels()
    test_vev_virt_different_labels()
    test_no_contraction_across_spaces()
    test_factor_preserved_through_normal_order()
