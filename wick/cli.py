#!/usr/bin/env python
"""Wick's Theorem Engine CLI."""

import sys
import argparse
from wick.helper import WickHelper
from wick.step_1_expand import expand_general_labels
from wick.step_2_normal_order import normal_order_fermi_vacuum
from wick.step_3_filter import filter_fully_contracted
from wick.step_4_deltas import apply_deltas
from wick.step_4b_reclassify import reclassify_occ_repulsion
from wick.step_5_labels import canonicalize_labels_pdaggerq_style
from wick.step_6_cleanup import cancel_terms
from wick.step_6d_graph_isomorphism import cancel_terms_graph_isomorphic
from wick.step_6a_linked_cluster import filter_unlinked_diagrams
from wick.step_6c_graph_connected import filter_connected_terms_graph
from wick.step_7_output import format_strings
from wick.step_8_denominator import extract_denominator
from wick.term import Term, Tensor


def print_section(title):
    """Print a section header."""
    print(f"\n{'=' * 72}")
    print(f" {title}")
    print(f"{'=' * 72}\n")


def print_terms(terms, limit=None):
    """Pretty-print a list of terms."""
    if not terms:
        print("  (empty)")
        return

    shown = terms if limit is None else terms[:limit]
    for i, term in enumerate(shown):
        factor_str = f"{term.factor:+.6f}" if term.sign > 0 else f"{-term.factor:+.6f}"
        ops_str = ", ".join(f"{op.label}{'†' if op.dagger else ''}" for op in term.operators)
        deltas_str = ", ".join(f"δ({d.i},{d.j})" for d in term.deltas)
        tensors_str = ", ".join(t.to_string() for t in term.tensors)

        print(f"  [{i}] factor={factor_str}, sign={term.sign:+d}")
        if tensors_str:
            print(f"       tensors: [{tensors_str}]")
        if ops_str:
            print(f"       ops: [{ops_str}]")
        if deltas_str:
            print(f"       deltas: [{deltas_str}]")

    if limit is not None and len(terms) > limit:
        print(f"  ... and {len(terms) - limit} more terms")
    print()


def print_strings(strings_output, limit=None):
    """Pretty-print strings() output format."""
    if not strings_output:
        print("  (empty)")
        return

    shown = strings_output if limit is None else strings_output[:limit]
    for i, row in enumerate(shown):
        factor = row[0]
        tensors = " * ".join(row[1:]) if len(row) > 1 else "(scalar)"
        print(f"  [{i:>2}] {factor:>8} * {tensors}")

    if limit is not None and len(strings_output) > limit:
        print(f"  ... and {len(strings_output) - limit} more terms")
    print()


def print_strings_with_denominators(strings_output, denominators, limit=None):
    """Pretty-print strings() output with energy denominators."""
    if not strings_output:
        print("  (empty)")
        return

    if limit is None:
        rows = strings_output
        dens = denominators
    else:
        rows = strings_output[:limit]
        dens = denominators[:limit]

    for i, (row, denom) in enumerate(zip(rows, dens)):
        factor = row[0]
        tensors = ' * '.join(row[1:]) if len(row) > 1 else "(scalar)"
        print(f"  [{i:>2}] {factor:>8} * {tensors}")
        print(f"         / ({denom})")

    if limit is not None and len(strings_output) > limit:
        print(f"  ... and {len(strings_output) - limit} more terms")
    print()


def pdaggerq_reference_strings(num_v):
    """Return pdaggerq strings() for an operator product of num_v 'v' terms."""
    import pdaggerq

    pq = pdaggerq.pq_helper("fermi")
    pq.set_left_operators([['1']])
    pq.add_operator_product(1.0, ['v'] * num_v)
    pq.simplify()
    return pq.strings()


def run_pipeline_verbose(wh, graph_connected_only=False, use_graph_simplification=True):
    """Run full pipeline with step-by-step verbose output."""
    print_section("Step 1: General Label Expansion")
    terms = list(wh._input_terms)
    counter = [0]
    terms = expand_general_labels(terms, counter)
    print(f"Output: {len(terms)} terms")
    print_terms(terms, limit=3)

    print_section("Step 2: Normal Ordering (Wick's Theorem)")
    terms = normal_order_fermi_vacuum(terms)
    print(f"Output: {len(terms)} terms")
    print_terms(terms, limit=3)

    print_section("Step 3: Filter Fully Contracted")
    terms = filter_fully_contracted(terms)
    print(f"Output: {len(terms)} terms")
    print_terms(terms, limit=3)

    print_section("Step 4: Delta Elimination")
    terms = apply_deltas(terms)
    print(f"Output: {len(terms)} terms")
    print_terms(terms, limit=3)

    print_section("Step 4b: Reclassify occ_repulsion")
    terms = reclassify_occ_repulsion(terms)
    print(f"Output: {len(terms)} terms")
    print_terms(terms, limit=3)

    print_section("Step 5: Label Canonicalization")
    terms = canonicalize_labels_pdaggerq_style(terms)
    print(f"Output: {len(terms)} terms")
    print_terms(terms, limit=3)

    print_section("Step 6: Consolidation")
    if use_graph_simplification:
        terms = cancel_terms_graph_isomorphic(terms)
    else:
        terms = cancel_terms(terms)
    print(f"Output: {len(terms)} terms")
    print_terms(terms, limit=3)

    if graph_connected_only:
        print_section("Step 6c: Graph Connected Filter")
        terms = filter_connected_terms_graph(terms)
        print(f"Output: {len(terms)} terms")
        print_terms(terms, limit=3)

    print_section("Step 6a: Linked-Diagram Filter")
    terms = filter_unlinked_diagrams(terms)
    print(f"Output: {len(terms)} terms")
    print_terms(terms, limit=3)

    print_section("Step 7: Format Strings")
    result = format_strings(terms)
    print(f"Output: {len(result)} terms")
    print_strings(result)

    print_section("Step 8: Energy Denominator Extraction")
    denoms = [extract_denominator(term) for term in terms]
    print(f"Output: {len(denoms)} denominators")
    print_strings_with_denominators(result, denoms)
    return result, denoms


def run_mp12(verbose=False, graph_connected_only=False, use_graph_simplification=True, show_reference=False):
    """Run MP1 + MP2 together in a single simplify() pass."""
    print_section("Running MP1+MP2 (Combined)")

    # Combined reference: strings from ['v'] and ['v','v'] together.
    ref_strings = None
    if show_reference:
        try:
            ref_strings = sorted(pdaggerq_reference_strings(1) + pdaggerq_reference_strings(2))
            print_section("Reference (pdaggerq)")
            print(f"pdaggerq combined output ({len(ref_strings)} terms):")
            print_strings(ref_strings)
        except Exception as exc:
            print_section("Reference (pdaggerq)")
            print(f"Could not compute pdaggerq reference: {exc}")

    wh = WickHelper(
        filter_unlinked=True,
        graph_connected_only=graph_connected_only,
        use_graph_simplification=use_graph_simplification,
    )
    wh.add_operator_product(1.0, ["v"])
    wh.add_operator_product(1.0, ["v", "v"])

    print_section("Running Pipeline")
    if verbose:
        result, denoms = run_pipeline_verbose(
            wh,
            graph_connected_only=graph_connected_only,
            use_graph_simplification=use_graph_simplification,
        )
        result = sorted(result)
    else:
        wh.simplify()
        result = sorted(wh.strings())
        denoms = [extract_denominator(t) for t in sorted(wh._result, key=lambda x: str(format_strings([x])[0]))]
    print(f"Our combined output ({len(result)} terms) with energy denominators:")
    print_strings_with_denominators(result, denoms)
    if ref_strings is not None:
        print(f"Exact string match vs pdaggerq: {result == ref_strings}")
    return result


def run_mp123(verbose=False, graph_connected_only=False, use_graph_simplification=True, show_reference=False):
    """Run MP1 + MP2 + MP3 together in a single simplify() pass."""
    print_section("Running MP1+MP2+MP3 (Combined)")

    # Combined reference from individual products.
    ref_strings = None
    if show_reference:
        try:
            ref_strings = sorted(
                pdaggerq_reference_strings(1)
                + pdaggerq_reference_strings(2)
                + pdaggerq_reference_strings(3)
            )
            print_section("Reference (pdaggerq)")
            print(f"pdaggerq combined output ({len(ref_strings)} terms):")
            print_strings(ref_strings)
        except Exception as exc:
            print_section("Reference (pdaggerq)")
            print(f"Could not compute pdaggerq reference: {exc}")

    wh = WickHelper(
        filter_unlinked=True,
        graph_connected_only=graph_connected_only,
        use_graph_simplification=use_graph_simplification,
    )
    wh.add_operator_product(1.0, ["v"])
    wh.add_operator_product(1.0, ["v", "v"])
    wh.add_operator_product(1.0, ["v", "v", "v"])

    print_section("Running Pipeline")
    if verbose:
        result, denoms = run_pipeline_verbose(
            wh,
            graph_connected_only=graph_connected_only,
            use_graph_simplification=use_graph_simplification,
        )
        result = sorted(result)
    else:
        wh.simplify()
        result = sorted(wh.strings())
        denoms = [extract_denominator(t) for t in sorted(wh._result, key=lambda x: str(format_strings([x])[0]))]
    print(f"Our combined output ({len(result)} terms) with energy denominators:")
    print_strings_with_denominators(result, denoms)
    if ref_strings is not None:
        print(f"Exact string match vs pdaggerq: {result == ref_strings}")
    return result


def run_example(
    name,
    factor,
    tensors,
    ops,
    verbose=False,
    use_v_product=False,
    graph_connected_only=False,
    use_graph_simplification=True,
    show_reference=False,
):
    """Run an example through the full pipeline with optional verbosity."""
    print_section(f"Running {name.upper()}")

    print(f"Input:")
    print(f"  factor = {factor}")
    print(f"  tensors = {tensors}")
    print(f"  ops = {ops}")

    # Show reference result first (for standard MPn examples)
    ref_strings = None
    if show_reference and name.lower().startswith("mp"):
        try:
            ref_strings = pdaggerq_reference_strings(len(tensors))
            print_section("Reference (pdaggerq)")
            print(f"pdaggerq output ({len(ref_strings)} terms):")
            print_strings(sorted(ref_strings))
        except Exception as exc:
            print_section("Reference (pdaggerq)")
            print(f"Could not compute pdaggerq reference: {exc}")

    # Create the term
    # Default CLI view: linked contributions with no energy-subspace projection.
    wh = WickHelper(
        filter_unlinked=True,
        graph_connected_only=graph_connected_only,
        use_graph_simplification=use_graph_simplification,
    )
    if use_v_product:
        wh.add_operator_product(1.0, ["v"] * len(tensors))
    else:
        wh.add_term(factor, tensors, ops)

    if not verbose:
        # Just run it
        print_section("Running Pipeline")
        wh.simplify()
        result = wh.strings()
        denoms = wh.denominators()
        print(f"Our output ({len(result)} terms) with energy denominators:")
        print_strings_with_denominators(result, denoms)
        if ref_strings is not None:
            print(f"Exact string match vs pdaggerq: {sorted(result) == sorted(ref_strings)}")
        return result

    result, _ = run_pipeline_verbose(
        wh,
        graph_connected_only=graph_connected_only,
        use_graph_simplification=use_graph_simplification,
    )
    if ref_strings is not None:
        print(f"Exact string match vs pdaggerq: {sorted(result) == sorted(ref_strings)}")

    return result


def run_mpn(order, verbose=False, graph_connected_only=False, use_graph_simplification=True, show_reference=False):
    """Run generic MPn using n copies of the fluctuation potential operator."""
    if order < 1:
        raise ValueError("MP order must be >= 1")

    print_section(f"Running MP{order}")
    print("Input:")
    print(f"  operator product = ['v'] * {order}")

    ref_strings = None
    if show_reference:
        try:
            ref_strings = sorted(pdaggerq_reference_strings(order))
            print_section("Reference (pdaggerq)")
            print(f"pdaggerq output ({len(ref_strings)} terms):")
            print_strings(ref_strings)
        except Exception as exc:
            print_section("Reference (pdaggerq)")
            print(f"Could not compute pdaggerq reference: {exc}")

    wh = WickHelper(
        filter_unlinked=True,
        graph_connected_only=graph_connected_only,
        use_graph_simplification=use_graph_simplification,
    )
    wh.add_operator_product(1.0, ["v"] * order)

    print_section("Running Pipeline")
    if verbose:
        result, denoms = run_pipeline_verbose(
            wh,
            graph_connected_only=graph_connected_only,
            use_graph_simplification=use_graph_simplification,
        )
    else:
        wh.simplify()
        result = wh.strings()
        denoms = wh.denominators()
        print(f"Our output ({len(result)} terms) with energy denominators:")
        print_strings_with_denominators(result, denoms)

    if ref_strings is not None:
        print(f"Exact string match vs pdaggerq: {sorted(result) == sorted(ref_strings)}")
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Wick's Theorem Engine CLI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  wick --mp1
  wick --mp2 --verbose
  wick --mp 5 --verbose
  wick --mp12
  wick --mp123 --reference
"""
    )

    parser.add_argument('--mp1', action='store_true', help='Run MP1 example')
    parser.add_argument('--mp2', action='store_true', help='Run MP2 example')
    parser.add_argument('--mp', type=int, default=None, help='Run generic MPn (e.g. --mp 5)')
    parser.add_argument('--mp12', action='store_true', help='Run MP1+MP2 combined example')
    parser.add_argument('--mp123', action='store_true', help='Run MP1+MP2+MP3 combined example')
    parser.add_argument('--mp3', action='store_true', help='Run MP3 example')
    parser.add_argument('--mp4', action='store_true', help='Run MP4 example (slow, 65K+ terms)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Show detailed output for each pipeline step')
    parser.add_argument('--graph-connected', action='store_true',
                       help='Enable graph-based connected-term filtering after cleanup')
    parser.add_argument('--graph-simplify', action='store_true',
                       help='Use graph-isomorphic simplification (experimental)')
    parser.add_argument('--reference', action='store_true',
                       help='Show pdaggerq reference output alongside our output')

    args = parser.parse_args()

    # Determine what to run
    if not any([args.mp1, args.mp2, args.mp12, args.mp123, args.mp3, args.mp4, args.mp is not None]):
        parser.print_help()
        return

    if args.mp is not None:
        run_mpn(
            args.mp,
            verbose=args.verbose,
            graph_connected_only=args.graph_connected,
            use_graph_simplification=args.graph_simplify,
            show_reference=args.reference,
        )
    elif args.mp1:
        run_example("MP1", 1/4, ['g(p,q,r,s)'], ['+p','+q','-r','-s'], args.verbose, use_v_product=True, graph_connected_only=args.graph_connected, use_graph_simplification=args.graph_simplify, show_reference=args.reference)
    elif args.mp2:
        run_example("MP2", 1/16, ['g(p,q,r,s)', 'g(t,u,v,w)'],
                   ['+p','+q','-r','-s','+t','+u','-v','-w'], args.verbose, use_v_product=True, graph_connected_only=args.graph_connected, use_graph_simplification=args.graph_simplify, show_reference=args.reference)
    elif args.mp12:
        run_mp12(verbose=args.verbose, graph_connected_only=args.graph_connected, use_graph_simplification=args.graph_simplify, show_reference=args.reference)
    elif args.mp123:
        run_mp123(verbose=args.verbose, graph_connected_only=args.graph_connected, use_graph_simplification=args.graph_simplify, show_reference=args.reference)
    elif args.mp3:
        run_example("MP3", 1/64,
                   ['g(p,q,r,s)', 'g(t,u,v,w)', 'g(x,y,z,o)'],
                   ['+p','+q','-r','-s', '+t','+u','-v','-w', '+x','+y','-z','-o'],
                   args.verbose, use_v_product=True, graph_connected_only=args.graph_connected, use_graph_simplification=args.graph_simplify, show_reference=args.reference)
    elif args.mp4:
        import time
        print_section("MP4 - WARNING: This is slow (65K+ initial terms)")
        print("Starting MP4 pipeline. This may take several minutes...\n")
        t_start = time.time()
        try:
            run_example("MP4", 1/256,
                       ['g(p,q,r,s)', 'g(t,u,v,w)', 'g(x,y,z,o)', 'g(p1,q1,r1,s1)'],
                       ['+p','+q','-r','-s', '+t','+u','-v','-w', '+x','+y','-z','-o', '+p1','+q1','-r1','-s1'],
                       args.verbose, use_v_product=True, graph_connected_only=args.graph_connected, use_graph_simplification=args.graph_simplify, show_reference=args.reference)
            elapsed = time.time() - t_start
            print_section("Timing")
            print(f"MP4 completed in {elapsed:.1f} seconds")
        except KeyboardInterrupt:
            print("\n\nInterrupted by user.")
            sys.exit(1)
if __name__ == '__main__':
    main()
