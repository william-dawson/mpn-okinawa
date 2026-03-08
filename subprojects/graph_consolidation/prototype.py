from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
import re
from typing import Dict, Iterable, List, Sequence, Tuple

ERI_RE = re.compile(r"^<([^,]+),([^|]+)\|\|([^,]+),([^>]+)>$")
_OCC = set("ijklmn")
_VIRT = set("abcdef")


def index_space(idx: str) -> str:
    if idx and idx[0] in _OCC:
        return "occ"
    if idx and idx[0] in _VIRT:
        return "virt"
    return "gen"


@dataclass(frozen=True)
class ERITerm:
    factor: float
    tensors: Tuple[Tuple[str, str, str, str], ...]


@dataclass
class GraphView:
    term: ERITerm
    labels: Dict[str, List[Tuple[int, int]]]


def parse_row(row: Sequence[str]) -> ERITerm:
    factor = float(row[0])
    tensors: List[Tuple[str, str, str, str]] = []
    for t in row[1:]:
        m = ERI_RE.match(t)
        if not m:
            raise ValueError(f"Bad ERI string: {t}")
        tensors.append((m.group(1), m.group(2), m.group(3), m.group(4)))
    return ERITerm(factor=factor, tensors=tuple(tensors))


def build_graph(term: ERITerm) -> GraphView:
    labels: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for ti, tensor in enumerate(term.tensors):
        for slot, lbl in enumerate(tensor):
            labels[lbl].append((ti, slot))
    return GraphView(term=term, labels=dict(labels))


def connected_components(term: ERITerm) -> List[List[int]]:
    gv = build_graph(term)
    n = len(term.tensors)
    if n == 0:
        return []

    adj = {i: set() for i in range(n)}
    for incidences in gv.labels.values():
        tensors = [ti for ti, _ in incidences]
        for i in tensors:
            for j in tensors:
                if i != j:
                    adj[i].add(j)

    seen = set()
    comps = []
    for i in range(n):
        if i in seen:
            continue
        stack = [i]
        comp = []
        seen.add(i)
        while stack:
            cur = stack.pop()
            comp.append(cur)
            for nxt in adj[cur]:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        comps.append(sorted(comp))
    return comps


def is_connected(term: ERITerm) -> bool:
    comps = connected_components(term)
    return len(comps) <= 1


def mp3_topology_class(term: ERITerm) -> str:
    labels = {x for tensor in term.tensors for x in tensor}
    n_occ = sum(1 for x in labels if index_space(x) == "occ")
    n_virt = sum(1 for x in labels if index_space(x) == "virt")

    if not is_connected(term):
        return "disconnected"
    if n_occ == 2 and n_virt == 4:
        return "pp"
    if n_occ == 4 and n_virt == 2:
        return "hh"
    if n_occ == 3 and n_virt == 3:
        return "ph"
    return "other"


def canonical_graph_key(term: ERITerm, max_iter: int = 8) -> str:
    """
    WL-style bipartite refinement key (prototype, not a formal isomorphism certificate).

    Node sets:
    - tensor nodes: T{idx}
    - label nodes: L:{label}
    """
    gv = build_graph(term)

    t_colors = {i: "T:g" for i in range(len(term.tensors))}
    l_colors = {lbl: f"L:{index_space(lbl)}" for lbl in gv.labels}

    for _ in range(max_iter):
        new_l = {}
        for lbl, incidences in gv.labels.items():
            neigh = sorted(f"{t_colors[ti]}@{slot}" for ti, slot in incidences)
            new_l[lbl] = f"{l_colors[lbl]}|{";".join(neigh)}"

        new_t = {}
        for ti, tensor in enumerate(term.tensors):
            neigh = [f"{new_l[tensor[slot]]}@{slot}" for slot in range(4)]
            new_t[ti] = f"{t_colors[ti]}|{";".join(sorted(neigh))}"

        if new_l == l_colors and new_t == t_colors:
            break
        l_colors, t_colors = new_l, new_t

    t_hist = Counter(t_colors.values())
    l_hist = Counter(l_colors.values())

    t_part = "|".join(f"{k}#{v}" for k, v in sorted(t_hist.items()))
    l_part = "|".join(f"{k}#{v}" for k, v in sorted(l_hist.items()))
    return f"T[{t_part}]::L[{l_part}]"


def group_rows(rows: Iterable[Sequence[str]]) -> Dict[str, Dict[str, float]]:
    """Group rows by topology and graph key with summed factors."""
    grouped: Dict[str, Dict[str, float]] = defaultdict(lambda: defaultdict(float))
    for row in rows:
        t = parse_row(row)
        topo = mp3_topology_class(t)
        key = canonical_graph_key(t)
        grouped[topo][key] += t.factor
    return grouped


def summarize_rows(rows: Iterable[Sequence[str]]) -> Dict[str, int]:
    counts = Counter()
    for row in rows:
        t = parse_row(row)
        counts[mp3_topology_class(t)] += 1
    return dict(counts)
