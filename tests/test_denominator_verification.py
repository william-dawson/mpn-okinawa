"""
Verification of denominator extraction against theoretical expectations.

MP2 correlation energy formula:
  E_2 = sum_{i<j,a<b} [(V^ab_ij)^2 - V^ab_ij V^ba_ij] / (e_i + e_j - e_a - e_b)

The Wick's theorem expansion produces multiple diagram types, each with different
intermediate excitation structures. Our denominator extraction should reflect these.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pdaggerq
from wick.helper import WickHelper


def test_mp2_denominators_vs_pdaggerq():
    """Compare our MP2 output with pdaggerq's."""
    print("\n" + "="*70)
    print("MP2: Our code vs pdaggerq")
    print("="*70)

    # Our implementation
    wh = WickHelper()
    wh.add_term(1/16, ['g(p,q,r,s)', 'g(t,u,v,w)'],
                ['+p','+q','-s','-r','+t','+u','-w','-v'])
    wh.simplify()

    our_strings = wh.strings()
    our_denoms = wh.denominators()
    our_pairs = wh.terms_with_denominators()

    # pdaggerq implementation
    pq = pdaggerq.pq_helper("fermi")
    pq.set_left_operators([['1']])
    pq.add_operator_product(1.0, ['v', 'v'])
    pq.simplify()
    pq_strings = pq.strings()

    print(f"\nOur output ({len(our_strings)} terms):")
    for i, (strings_row, denom) in enumerate(our_pairs):
        factor = strings_row[0]
        tensors = ' * '.join(strings_row[1:])
        print(f"  [{i}] {factor:>8} * {tensors:30} / ({denom})")

    print(f"\npdaggerq output ({len(pq_strings)} terms):")
    for i, row in enumerate(pq_strings):
        factor = row[0]
        tensors = ' * '.join(row[1:])
        print(f"  [{i}] {factor:>8} * {tensors:30}")

    print(f"\nTerm counts: Ours={len(our_strings)}, pdaggerq={len(pq_strings)}")
    print(f"Denominator structures found: {len(set(our_denoms))} unique")
    for denom in sorted(set(our_denoms)):
        count = our_denoms.count(denom)
        print(f"  • {denom:30} appears {count}x")

    # Verify string magnitudes match (can't compare exact order due to canonicalization)
    our_factors = sorted([float(row[0]) for row in our_strings])
    pq_factors = sorted([float(row[0]) for row in pq_strings])

    print(f"\nFactor comparison (sorted by magnitude):")
    print(f"  Ours:     {our_factors}")
    print(f"  pdaggerq: {pq_factors}")
    print(f"  Match: {abs(sum(our_factors) - sum(pq_factors)) < 1e-6}")


def test_mp3_denominators_vs_pdaggerq():
    """Compare our MP3 output with pdaggerq's."""
    print("\n" + "="*70)
    print("MP3: Our code vs pdaggerq")
    print("="*70)

    # Our implementation
    wh = WickHelper()
    wh.add_term(1/64,
                ['g(p,q,r,s)', 'g(t,u,v,w)', 'g(A,B,C,D)'],
                ['+p','+q','-s','-r', '+t','+u','-w','-v', '+A','+B','-D','-C'])
    wh.simplify()

    our_strings = wh.strings()
    our_denoms = wh.denominators()
    our_pairs = wh.terms_with_denominators()

    # pdaggerq implementation
    pq = pdaggerq.pq_helper("fermi")
    pq.set_left_operators([['1']])
    pq.add_operator_product(1.0, ['v', 'v', 'v'])
    pq.simplify()
    pq_strings = pq.strings()

    print(f"\nOur output ({len(our_strings)} terms, showing first 5):")
    for i, (strings_row, denom) in enumerate(our_pairs[:5]):
        factor = strings_row[0]
        tensors = ' * '.join(strings_row[1:])
        print(f"  [{i}] {factor:>8} * {tensors:40} / ({denom})")
    if len(our_pairs) > 5:
        print(f"  ... and {len(our_pairs) - 5} more terms")

    print(f"\npdaggerq output ({len(pq_strings)} terms, showing first 5):")
    for i, row in enumerate(pq_strings[:5]):
        factor = row[0]
        tensors = ' * '.join(row[1:])
        print(f"  [{i}] {factor:>8} * {tensors:40}")
    if len(pq_strings) > 5:
        print(f"  ... and {len(pq_strings) - 5} more terms")

    print(f"\nTerm counts: Ours={len(our_strings)}, pdaggerq={len(pq_strings)}")
    print(f"Denominator structures found: {len(set(our_denoms))} unique")

    # Show distribution
    denom_counts = {}
    for d in our_denoms:
        denom_counts[d] = denom_counts.get(d, 0) + 1

    print("Denominator distribution:")
    for denom in sorted(denom_counts.keys()):
        count = denom_counts[denom]
        print(f"  • {denom:40} : {count}x")

    # Verify factors
    our_factors = sorted([float(row[0]) for row in our_strings])
    pq_factors = sorted([float(row[0]) for row in pq_strings])

    print(f"\nFactor sum comparison:")
    print(f"  Ours:     {sum(our_factors):.6f}")
    print(f"  pdaggerq: {sum(pq_factors):.6f}")
    print(f"  Match: {abs(sum(our_factors) - sum(pq_factors)) < 1e-6}")


if __name__ == '__main__':
    test_mp2_denominators_vs_pdaggerq()
    test_mp3_denominators_vs_pdaggerq()
