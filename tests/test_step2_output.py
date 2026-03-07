"""
Step 2 isolated test: output formatting.

Hand-craft Terms matching pdaggerq's MP2 output and verify format_strings()
produces exactly the pdaggerq string format.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from wick.term import Term, Tensor
from wick.step_7_output import format_strings, _fmt_factor


def test_fmt_factor():
    # Expected values verified against actual pdaggerq output
    cases = [
        ( 0.25,   '+0.250'),
        (-0.25,   '-0.250'),
        ( 0.125,  '+0.1250'),
        ( 0.375,  '+0.3750'),
        ( 1.0,    '+1.00'),
        (-1.5,    '-1.50'),
    ]
    for value, expected in cases:
        got = _fmt_factor(value)
        assert got == expected, f"_fmt_factor({value}): got {got!r}, want {expected!r}"
    print("PASS: _fmt_factor")


def test_format_strings_mp2():
    """Match the exact two terms pdaggerq produces for MP2."""
    # ['+0.250', '<k,i||k,i>', '<l,j||l,j>']
    t1 = Term(0.25, +1,
              [Tensor('g', ['k','i','k','i']),
               Tensor('g', ['l','j','l','j'])],
              operators=[])

    # ['+0.250', '<a,b||j,i>', '<j,i||a,b>']
    t2 = Term(0.25, +1,
              [Tensor('g', ['a','b','j','i']),
               Tensor('g', ['j','i','a','b'])],
              operators=[])

    result = format_strings([t1, t2])

    assert result[0] == ['+0.250', '<k,i||k,i>', '<l,j||l,j>'], result[0]
    assert result[1] == ['+0.250', '<a,b||j,i>', '<j,i||a,b>'], result[1]
    print("PASS: format_strings (MP2 terms)")


def test_negative_factor():
    t = Term(0.125, -1, [Tensor('g', ['i','j','a','b'])], operators=[])
    result = format_strings([t])
    assert result[0][0] == '-0.1250', result[0][0]  # 0.125 -> 4 decimal places
    print("PASS: format_strings (negative factor)")


if __name__ == '__main__':
    test_fmt_factor()
    test_format_strings_mp2()
    test_negative_factor()
