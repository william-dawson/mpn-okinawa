"""
Step 5 isolated test: contraction filter.

filter_fully_contracted() keeps only terms with no remaining operators.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from wick.term import Term, Tensor, Operator
from wick.step_3_filter import filter_fully_contracted


def test_keeps_fully_contracted():
    t = Term(0.25, 1, [Tensor('g', ['i','j','a','b'])], operators=[])
    result = filter_fully_contracted([t])
    assert result == [t]
    print("PASS: fully contracted term kept")


def test_drops_term_with_operators():
    op = Operator('i', True)
    t = Term(0.25, 1, [Tensor('g', ['i','j','a','b'])], operators=[op])
    result = filter_fully_contracted([t])
    assert result == []
    print("PASS: term with operators dropped")


def test_mixed_list():
    t_contracted   = Term(0.25, 1, [Tensor('g', ['i','j','a','b'])], operators=[])
    t_uncontracted = Term(0.25, 1, [Tensor('g', ['i','j','a','b'])], operators=[Operator('i', False)])
    result = filter_fully_contracted([t_contracted, t_uncontracted])
    assert result == [t_contracted]
    assert len(result) == 1
    print("PASS: mixed list — only contracted term kept")


def test_empty_list():
    assert filter_fully_contracted([]) == []
    print("PASS: empty list")


def test_all_uncontracted():
    terms = [
        Term(1.0, 1, [], operators=[Operator('a', True)]),
        Term(1.0, 1, [], operators=[Operator('i', False)]),
    ]
    assert filter_fully_contracted(terms) == []
    print("PASS: all uncontracted — empty result")


if __name__ == '__main__':
    test_keeps_fully_contracted()
    test_drops_term_with_operators()
    test_mixed_list()
    test_empty_list()
    test_all_uncontracted()
