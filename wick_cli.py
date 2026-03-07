#!/usr/bin/env python
"""
Wick's Theorem Engine CLI

Interactive command-line interface for testing and exploring the pipeline.

Usage:
    python wick_cli.py                          # Start interactive mode
    python wick_cli.py --mp1                    # Run MP1 example
    python wick_cli.py --mp2                    # Run MP2 example
    python wick_cli.py --mp3                    # Run MP3 example
    python wick_cli.py --verbose --mp2          # Run with detailed output at each step
"""

import sys
import argparse
from wick.helper import WickHelper
from wick.step_1_expand import expand_general_labels
from wick.step_2_normal_order import normal_order_fermi_vacuum
from wick.step_3_filter import filter_fully_contracted
from wick.step_4_deltas import apply_deltas
from wick.step_5_labels import canonicalize_labels
from wick.step_6_cleanup import cancel_terms
from wick.step_7_output import format_strings
from wick.step_8_denominator import extract_denominator
from wick.term import Term, Tensor


def print_section(title):
    """Print a section header."""
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}\n")


def print_terms(terms, limit=10):
    """Pretty-print a list of terms."""
    if not terms:
        print("  (empty)")
        return

    for i, term in enumerate(terms[:limit]):
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

    if len(terms) > limit:
        print(f"  ... and {len(terms) - limit} more terms")
    print()


def print_strings(strings_output, limit=10):
    """Pretty-print strings() output format."""
    if not strings_output:
        print("  (empty)")
        return

    for i, row in enumerate(strings_output[:limit]):
        print(f"  [{i}] {row}")

    if len(strings_output) > limit:
        print(f"  ... and {len(strings_output) - limit} more terms")
    print()


def print_strings_with_denominators(strings_output, denominators, limit=10):
    """Pretty-print strings() output with energy denominators."""
    if not strings_output:
        print("  (empty)")
        return

    for i, (row, denom) in enumerate(zip(strings_output[:limit], denominators[:limit])):
        factor = row[0]
        tensors = ' * '.join(row[1:]) if len(row) > 1 else "(scalar)"
        print(f"  [{i}] {factor:>8} * {tensors}")
        print(f"         / ({denom})")

    if len(strings_output) > limit:
        print(f"  ... and {len(strings_output) - limit} more terms")
    print()


def run_example(name, factor, tensors, ops, verbose=False):
    """Run an example through the full pipeline with optional verbosity."""
    print_section(f"Running {name.upper()}")

    print(f"Input:")
    print(f"  factor = {factor}")
    print(f"  tensors = {tensors}")
    print(f"  ops = {ops}")

    # Create the term
    # Keep CLI examples aligned with pdaggerq strings() reference behavior.
    wh = WickHelper(filter_unlinked=False)
    wh.add_term(factor, tensors, ops)

    if not verbose:
        # Just run it
        print_section("Running Pipeline")
        wh.simplify()
        result = wh.strings()
        denoms = wh.denominators()
        print(f"Output ({len(result)} terms) with energy denominators:")
        print_strings_with_denominators(result, denoms)
        return result

    # Verbose mode: step through each stage
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

    print_section("Step 5: Label Canonicalization")
    terms = canonicalize_labels(terms)
    print(f"Output: {len(terms)} terms")
    print_terms(terms, limit=3)

    print_section("Step 6: Cleanup & ERI Canonicalization")
    terms = cancel_terms(terms)
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

    return result


def interactive_mode():
    """Interactive command loop."""
    print_section("Wick's Theorem Engine - Interactive Mode")
    print("""Commands:
  mp1                 Run MP1 example
  mp2                 Run MP2 example
  mp3                 Run MP3 example
  mp4                 Run MP4 example (SLOW - 65K+ terms)
  custom              Enter custom input
  help                Show this help
  quit/exit           Exit
""")

    while True:
        try:
            cmd = input("wick> ").strip().lower()

            if not cmd:
                continue
            elif cmd == "mp1":
                run_example("MP1", 1/4, ['g(p,q,r,s)'], ['+p','+q','-r','-s'])
            elif cmd == "mp2":
                run_example("MP2", 1/16, ['g(p,q,r,s)', 'g(t,u,v,w)'],
                           ['+p','+q','-r','-s','+t','+u','-v','-w'])
            elif cmd == "mp3":
                run_example("MP3", 1/64,
                           ['g(p,q,r,s)', 'g(t,u,v,w)', 'g(x,y,z,o)'],
                           ['+p','+q','-r','-s', '+t','+u','-v','-w', '+x','+y','-z','-o'])
            elif cmd == "mp4":
                print("\nWARNING: MP4 is very slow (65K+ initial terms).")
                confirm = input("Continue? (y/n): ").strip().lower()
                if confirm == 'y':
                    import time
                    t_start = time.time()
                    run_example("MP4", 1/256,
                               ['g(p,q,r,s)', 'g(t,u,v,w)', 'g(x,y,z,o)', 'g(p1,q1,r1,s1)'],
                               ['+p','+q','-r','-s', '+t','+u','-v','-w', '+x','+y','-z','-o', '+p1','+q1','-r1','-s1'])
                    elapsed = time.time() - t_start
                    print_section("Timing")
                    print(f"MP4 completed in {elapsed:.1f} seconds")
            elif cmd == "custom":
                print("\nEnter custom input:")
                factor = float(input("  factor: "))
                tensors_str = input("  tensors (comma-separated, e.g. 'g(p,q,r,s), g(a,b,c,d)'): ")
                tensors = [t.strip() for t in tensors_str.split(",")]
                ops_str = input("  ops (comma-separated, e.g. '+p, +q, -s, -r'): ")
                ops = [o.strip() for o in ops_str.split(",")]

                run_example("Custom", factor, tensors, ops)
            elif cmd in ["help", "h", "?"]:
                print("""Commands:
  mp1                 Run MP1 example
  mp2                 Run MP2 example
  mp3                 Run MP3 example
  custom              Enter custom input
  help                Show this help
  quit/exit           Exit
""")
            elif cmd in ["quit", "exit", "q"]:
                print("Exiting.")
                sys.exit(0)
            else:
                print(f"Unknown command: {cmd}. Type 'help' for options.")
        except KeyboardInterrupt:
            print("\nExiting.")
            sys.exit(0)
        except Exception as e:
            print(f"Error: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Wick's Theorem Engine CLI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  python wick_cli.py --mp1                    # Run MP1
  python wick_cli.py --mp2 --verbose          # Run MP2 with step-by-step output
  python wick_cli.py --mp4                    # Run MP4 (slow!)
  python wick_cli.py                          # Interactive mode
"""
    )

    parser.add_argument('--mp1', action='store_true', help='Run MP1 example')
    parser.add_argument('--mp2', action='store_true', help='Run MP2 example')
    parser.add_argument('--mp3', action='store_true', help='Run MP3 example')
    parser.add_argument('--mp4', action='store_true', help='Run MP4 example (slow, 65K+ terms)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Show detailed output for each pipeline step')

    args = parser.parse_args()

    # Determine what to run
    if args.mp1:
        run_example("MP1", 1/4, ['g(p,q,r,s)'], ['+p','+q','-r','-s'], args.verbose)
    elif args.mp2:
        run_example("MP2", 1/16, ['g(p,q,r,s)', 'g(t,u,v,w)'],
                   ['+p','+q','-r','-s','+t','+u','-v','-w'], args.verbose)
    elif args.mp3:
        run_example("MP3", 1/64,
                   ['g(p,q,r,s)', 'g(t,u,v,w)', 'g(x,y,z,o)'],
                   ['+p','+q','-r','-s', '+t','+u','-v','-w', '+x','+y','-z','-o'],
                   args.verbose)
    elif args.mp4:
        import time
        print_section("MP4 - WARNING: This is slow (65K+ initial terms)")
        print("Starting MP4 pipeline. This may take several minutes...\n")
        t_start = time.time()
        try:
            run_example("MP4", 1/256,
                       ['g(p,q,r,s)', 'g(t,u,v,w)', 'g(x,y,z,o)', 'g(p1,q1,r1,s1)'],
                       ['+p','+q','-r','-s', '+t','+u','-v','-w', '+x','+y','-z','-o', '+p1','+q1','-r1','-s1'],
                       args.verbose)
            elapsed = time.time() - t_start
            print_section("Timing")
            print(f"MP4 completed in {elapsed:.1f} seconds")
        except KeyboardInterrupt:
            print("\n\nInterrupted by user.")
            sys.exit(1)
    else:
        # Interactive mode
        interactive_mode()


if __name__ == '__main__':
    main()
