"""
Step 6 isolated test: general label expansion.

expand_general_labels() replaces each general label (p, q, r, s, ...)
with either a fresh occupied label (o<n>) or a fresh virtual label (v<n>),
producing 2^n terms for n general labels.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from wick.term import Term, Tensor, Operator, is_occ, is_virt, is_general
from wick.step_1_expand import expand_general_labels


def _all_labels(term):
    """Collect all index labels from tensors and operators in a term."""
    labels = []
    for t in term.tensors:
        labels += t.indices
    for op in term.operators:
        labels.append(op.label)
    return labels


def test_no_general_labels_unchanged():
    t = Term(1.0, 1, [Tensor('g', ['i','j','a','b'])], [])
    result = expand_general_labels([t])
    assert len(result) == 1
    assert result[0].tensors[0].indices == ['i','j','a','b']
    print("PASS: no general labels — term unchanged")


def test_single_general_label_gives_two_terms():
    # g(p, i, a, b) with operator +p
    t = Term(1.0, 1,
             [Tensor('g', ['p','i','a','b'])],
             [Operator('p', True)])
    result = expand_general_labels([t])
    assert len(result) == 2
    # one term should have p -> occ, one -> virt
    spaces = set()
    for term in result:
        lbl = term.tensors[0].indices[0]
        assert not is_general(lbl), f"label {lbl!r} still general"
        # operator label must match tensor label
        assert term.operators[0].label == lbl
        spaces.add('occ' if is_occ(lbl) else 'virt')
    assert spaces == {'occ', 'virt'}
    print("PASS: single general label -> 2 terms with consistent labels")


def test_two_general_labels_gives_four_terms():
    t = Term(1.0, 1, [Tensor('g', ['p','q','i','j'])], [])
    result = expand_general_labels([t])
    assert len(result) == 4
    combos = set()
    for term in result:
        p_lbl = term.tensors[0].indices[0]
        q_lbl = term.tensors[0].indices[1]
        assert not is_general(p_lbl)
        assert not is_general(q_lbl)
        combos.add(('occ' if is_occ(p_lbl) else 'virt',
                    'occ' if is_occ(q_lbl) else 'virt'))
    assert combos == {('occ','occ'), ('occ','virt'), ('virt','occ'), ('virt','virt')}
    print("PASS: two general labels -> 4 terms covering all oo/ov/vo/vv combos")


def test_four_general_labels_gives_sixteen_terms():
    t = Term(1.0, 1, [Tensor('g', ['p','q','r','s'])],
             [Operator('p', True), Operator('q', True),
              Operator('s', False), Operator('r', False)])
    result = expand_general_labels([t])
    assert len(result) == 16
    print("PASS: four general labels -> 16 terms")


def test_fresh_labels_unique_within_term():
    # All four labels in a single term must map to distinct fresh labels
    t = Term(1.0, 1, [Tensor('g', ['p','q','r','s'])], [])
    result = expand_general_labels([t])
    for term in result:
        lbls = term.tensors[0].indices
        assert len(set(lbls)) == 4, f"duplicate fresh labels in {lbls}"
    print("PASS: fresh labels are distinct within each expanded term")


def test_fresh_labels_unique_across_terms():
    # Labels must be globally unique so terms don't share summation indices
    t = Term(1.0, 1, [Tensor('g', ['p','q','r','s'])], [])
    result = expand_general_labels([t])
    all_label_sets = [set(term.tensors[0].indices) for term in result]
    for i, s1 in enumerate(all_label_sets):
        for j, s2 in enumerate(all_label_sets):
            if i != j:
                assert s1.isdisjoint(s2), \
                    f"terms {i} and {j} share labels: {s1 & s2}"
    print("PASS: fresh labels are unique across all expanded terms")


def test_factor_sign_preserved():
    t = Term(0.0625, -1, [Tensor('g', ['p','q','r','s'])], [])
    result = expand_general_labels([t])
    for term in result:
        assert term.factor == 0.0625
        assert term.sign == -1
    print("PASS: factor and sign preserved across expansion")


def test_non_general_labels_untouched():
    # i,j,a,b are already specific — they must not be changed
    t = Term(1.0, 1, [Tensor('g', ['p','i','a','b'])], [])
    result = expand_general_labels([t])
    for term in result:
        indices = term.tensors[0].indices
        assert indices[1] == 'i'
        assert indices[2] == 'a'
        assert indices[3] == 'b'
    print("PASS: non-general labels untouched during expansion")


if __name__ == '__main__':
    test_no_general_labels_unchanged()
    test_single_general_label_gives_two_terms()
    test_two_general_labels_gives_four_terms()
    test_four_general_labels_gives_sixteen_terms()
    test_fresh_labels_unique_within_term()
    test_fresh_labels_unique_across_terms()
    test_factor_sign_preserved()
    test_non_general_labels_untouched()
