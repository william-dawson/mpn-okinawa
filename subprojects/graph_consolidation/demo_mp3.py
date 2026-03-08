from __future__ import annotations

import os
import sys
import pdaggerq

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from wick.helper import WickHelper
from subprojects.graph_consolidation.prototype import group_rows, summarize_rows


def pdaggerq_mp3_rows():
    pq = pdaggerq.pq_helper("fermi")
    pq.set_left_operators([["1"]])
    pq.add_operator_product(1.0, ["v", "v", "v"])
    pq.simplify()
    return sorted(pq.strings())


def ours_mp3_rows():
    wh = WickHelper(filter_unlinked=False, project_energy_subspace=False)
    wh.add_operator_product(1.0, ["v", "v", "v"])
    wh.simplify()
    return sorted(wh.strings())


def print_summary(name, rows):
    print(f"\n{name}")
    print(f"rows: {len(rows)}")
    print(f"topology counts: {summarize_rows(rows)}")

    grouped = group_rows(rows)
    for topo in sorted(grouped):
        coeffs = list(grouped[topo].values())
        print(f"  {topo}: groups={len(grouped[topo])}, coeff_sum={sum(coeffs):+.6f}")


def main():
    ref = pdaggerq_mp3_rows()
    ours = ours_mp3_rows()

    print_summary("pdaggerq MP3", ref)
    print_summary("ours MP3", ours)


if __name__ == "__main__":
    main()
