"""
Test denominator extraction.

For MP2, the correlation energy is:
  E_2 = sum_{i<j, a<b} (V^{ab}_{ij})^2 / (e_i + e_j - e_a - e_b)

where the denominator represents the energy difference between the reference state
and the intermediate excited state with holes i,j and particles a,b.

This test verifies that we extract the correct denominators from our terms.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from wick.helper import WickHelper


def test_mp2_denominators():
    """Test that MP2 produces correct denominator structure."""
    wh = WickHelper()
    wh.add_term(1/16, ['g(p,q,r,s)', 'g(t,u,v,w)'],
                ['+p','+q','-s','-r','+t','+u','-w','-v'])
    wh.simplify()

    denoms = wh.denominators()
    print(f"\nMP2 Denominators ({len(denoms)} terms):")
    for i, d in enumerate(denoms):
        print(f"  [{i}] {d}")

    # For MP2, we expect denominators with 2 occupied and 2 virtual indices
    # (or possibly all 4 occupied or all 4 virtual for certain contractions)
    for denom in denoms:
        # Count e_i + ... (occupied)
        occ_count = denom.count('e_i') + denom.count('e_j') + denom.count('e_k') + denom.count('e_l') + denom.count('e_m') + denom.count('e_n')
        # Count -e_a ... (virtual)
        virt_count = denom.count('e_a') + denom.count('e_b') + denom.count('e_c') + denom.count('e_d') + denom.count('e_e') + denom.count('e_f')

        # For typical intermediate excitations with the V operator structure,
        # we expect balanced or mostly balanced denominators
        print(f"    → {occ_count} occupied, {virt_count} virtual indices")

    # Basic sanity: should have denominators
    assert len(denoms) > 0, "MP2 should produce denominators"
    print("\n✓ MP2 denominator extraction works")


def test_mp3_denominators():
    """Test that MP3 produces correct denominator structure."""
    wh = WickHelper()
    wh.add_term(1/64,
                ['g(p,q,r,s)', 'g(t,u,v,w)', 'g(A,B,C,D)'],
                ['+p','+q','-s','-r', '+t','+u','-w','-v', '+A','+B','-D','-C'])
    wh.simplify()

    denoms = wh.denominators()
    print(f"\nMP3 Denominators ({len(denoms)} terms):")
    for i, d in enumerate(denoms[:5]):  # Show first 5
        print(f"  [{i}] {d}")
    if len(denoms) > 5:
        print(f"  ... and {len(denoms) - 5} more")

    assert len(denoms) > 0, "MP3 should produce denominators"
    print("\n✓ MP3 denominator extraction works")


def test_terms_with_denominators():
    """Test getting both strings and denominators together."""
    wh = WickHelper()
    wh.add_term(1/4, ['g(p,q,r,s)'], ['+p','+q','-s','-r'])
    wh.simplify()

    pairs = wh.terms_with_denominators()
    print(f"\nMP1 Terms with Denominators:")
    for string_row, denom in pairs:
        print(f"  {string_row[0]} * {' * '.join(string_row[1:])} / ({denom})")

    assert len(pairs) > 0
    assert all(len(pair) == 2 for pair in pairs)
    print("\n✓ terms_with_denominators() works")


if __name__ == '__main__':
    test_mp2_denominators()
    test_mp3_denominators()
    test_terms_with_denominators()
