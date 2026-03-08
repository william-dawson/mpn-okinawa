import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from wick.term import Term, Tensor, Operator, Delta, is_general, is_occ
from wick.step_1_expand import expand_general_labels
from wick.step_2_normal_order import normal_order_fermi_vacuum
from wick.step_3_filter import filter_fully_contracted
from wick.step_4_deltas import apply_deltas
from wick.step_4b_reclassify import reclassify_occ_repulsion
from wick.step_5_labels import canonicalize_labels_pdaggerq_style
from wick.step_6_cleanup import cancel_terms
from wick.step_6a_linked_cluster import filter_unlinked_diagrams


def test_expand_general_labels_basic():
    t = Term(1.0, 1, [Tensor('g', ['p', 'q', 'i', 'a'])], [Operator('p', True), Operator('q', False)])
    out = expand_general_labels([t])

    assert len(out) == 4
    for term in out:
        for tensor in term.tensors:
            for idx in tensor.indices:
                assert not is_general(idx)
        for op in term.operators:
            assert not is_general(op.label)


def test_normal_order_then_filter_fully_contracted():
    t = Term(1.0, 1, [Tensor('g', ['i', 'j', 'i', 'j'])], [Operator('i', True), Operator('j', False)])
    ordered = normal_order_fermi_vacuum([t])
    contracted = filter_fully_contracted(ordered)

    assert len(ordered) == 2
    assert len(contracted) == 1
    assert contracted[0].operators == []
    assert len(contracted[0].deltas) == 1


def test_apply_deltas_eliminates_delta_and_unifies_label_pair():
    t = Term(1.0, 1, [Tensor('g', ['o1', 'v1', 'o2', 'v2'])], [Operator('o1', True)], [Delta('o1', 'o2')])
    out = apply_deltas([t])[0]

    assert out.deltas == []
    labels_after = {idx for tensor in out.tensors for idx in tensor.indices if idx in {'o1', 'o2'}}
    assert len(labels_after) == 1
    assert out.operators[0].label in {'o1', 'o2'}


def test_reclassify_occ_repulsion_into_eri():
    t = Term(1.0, 1, [Tensor('occ_repulsion', ['o1', 'o2'])], [])
    out = reclassify_occ_repulsion([t])[0]

    assert len(out.tensors) == 1
    tensor = out.tensors[0]
    assert tensor.name == 'g'
    assert tensor.indices[0] == 'o1'
    assert tensor.indices[2] == 'o2'
    assert tensor.indices[1] == tensor.indices[3]
    assert is_occ(tensor.indices[1])


def test_canonicalize_labels_removes_fresh_prefixes():
    t = Term(1.0, 1, [Tensor('g', ['o3', 'v9', 'o4', 'v10'])], [])
    out = canonicalize_labels_pdaggerq_style([t])[0]

    for idx in out.tensors[0].indices:
        assert not (idx.startswith('o') and idx[1:].isdigit())
        assert not (idx.startswith('v') and idx[1:].isdigit())


def test_cancel_terms_sums_and_cancels_identical_terms():
    t1 = Term(0.25, 1, [Tensor('g', ['i', 'j', 'a', 'b'])], [])
    t2 = Term(0.25, -1, [Tensor('g', ['i', 'j', 'a', 'b'])], [])
    out = cancel_terms([t1, t2])
    assert out == []


def test_linked_filter_behaves_reasonably():
    connected = Term(1.0, 1, [Tensor('g', ['i', 'j', 'a', 'b']), Tensor('g', ['a', 'b', 'k', 'l'])], [])
    disconnected = Term(1.0, 1, [Tensor('g', ['i', 'j', 'i', 'j']), Tensor('g', ['k', 'l', 'k', 'l'])], [])
    linked = filter_unlinked_diagrams([connected, disconnected])
    assert linked == [connected]
