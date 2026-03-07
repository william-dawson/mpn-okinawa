# Project Plan: Wick's Theorem Engine

## Goal

Build a pure Python implementation of the second-quantized operator algebra
pipeline that pdaggerq implements in C++. Our code takes explicit operator
strings as input (rather than pdaggerq's high-level BCH/commutator API),
applies Wick's theorem with respect to the Fermi vacuum, and produces output
that matches pdaggerq's `strings()` format exactly.

**Target**: MP2 and MP3 numerators (fully contracted vacuum-to-vacuum matrix
elements of products of the antisymmetrized two-electron operator V).

---

## Non-Negotiable Rules

1. **Never advance without full confidence a step is done.**
2. **Every step is tested in isolation before being wired into the pipeline.**
3. **Every integration swap is confirmed passing before moving to the next step.**
4. **When confused about what pdaggerq does internally, READ THE SOURCE.**
   The reference implementation lives at `pdaggerq/pdaggerq/` (headers) and
   `pdaggerq/pdaggerq/*.py`. Do not guess. Read it.
5. **After each successful swap, stop and report to the user. Wait for feedback
   before proceeding.**

---

## Module Structure (Refactored by Pipeline Step)

### Data Structures
| Module | Purpose |
|--------|---------|
| `wick/term.py` | Core data classes: `Term`, `Tensor`, `Operator`, `Delta` |
| `wick/helper.py` | Main orchestrator: `WickHelper` class and pipeline runner |
| `wick/parse_pdaggerq.py` | Utilities for parsing pdaggerq output |

### Pipeline Steps (in execution order)
| Step | Module | Function | Input → Output |
|------|--------|----------|-----------------|
| 1 | `wick/step_1_expand.py` | `expand_general_labels()` | General indices (p,q,r,s) → occ/virt cases |
| 2 | `wick/step_2_normal_order.py` | `normal_order_fermi_vacuum()` | Operators → Fermi-vacuum normal order + deltas |
| 3 | `wick/step_3_filter.py` | `filter_fully_contracted()` | All terms → Fully-contracted terms only |
| 4 | `wick/step_4_deltas.py` | `apply_deltas()` | Deltas → Label substitutions |
| 5 | `wick/step_5_labels.py` | `canonicalize_labels()` | Fresh labels (o1,v2,...) → Conventional (i,j,.../a,b,...) |
| 6 | `wick/step_6_cleanup.py` | `cancel_terms()` | Terms with ERI canonicalization → Summed equivalent terms |
| 7 | `wick/step_7_output.py` | `format_strings()` | Terms → pdaggerq string format |

## Reference Code Locations (pdaggerq)

| What | Where |
|------|-------|
| pq_helper API | `pdaggerq/pdaggerq/pq_helper.h` |
| pq_string data structure | `pdaggerq/pdaggerq/pq_string.h` |
| Tensor types | `pdaggerq/pdaggerq/pq_tensor.h` |
| Normal ordering (Fermi vacuum) | `pdaggerq/pdaggerq/pq_swap_operators.h` |
| Delta/label utilities | `pdaggerq/pdaggerq/pq_utils.h` |
| Factor formatting | `pdaggerq/pdaggerq/pq_string.h` → `minimum_precision()` |
| Python examples | `pdaggerq/examples/ccsd_energy.py`, `ccsd_singles.py` |

---

## Input Convention

```python
wh.add_term(
    factor,                          # float coefficient
    tensors=['g(p,q,r,s)', ...],     # tensor strings
    ops=['+p', '+q', '-s', '-r'],    # operator strings
)
```

- `+label` = a†_label (true vacuum creator)
- `-label` = a_label  (true vacuum annihilator)
- `i,j,k,l,m,n` = occupied;  `a,b,c,d,e,f` = virtual;  anything else = general
- General labels are expanded into occ + virt cases before Wick's theorem

For MP2: `add_term(1/16, ['g(p,q,r,s)', 'g(t,u,v,w)'], ['+p','+q','-s','-r','+t','+u','-w','-v'])`
For MP3: same pattern with three V's, factor 1/64

---

## Pipeline

```
add_term(s)
    │
    ▼
[1] expand_general_labels      p,q,r,s → occ/virt cases  (2^n terms per general label)
    │
    ▼
[2] normal_order_fermi_vacuum  Wick's theorem: swap until normal ordered, insert deltas
    │
    ▼
[3] filter_fully_contracted    drop any term with remaining operators
    │
    ▼
[4] apply_deltas               substitute δ_{ij}: replace one label everywhere, remove delta
    │
    ▼
[5] canonicalize_labels        o1→i, o2→j, ..., v1→a, v2→b, ...
    │
    ▼
[6] cancel_terms               canonicalize integrals <p,q||r,s>, sum identical terms, drop zeros
    │
    ▼
[7] format_strings             → [[factor_str, tensor_str, ...], ...]
```

**New Step 6 Enhancement**: ERI Antisymmetry Canonicalization
- Sort bra pair (p,q): if p > q, swap and apply sign (-1)
- Sort ket pair (r,s): if r > s, swap and apply sign (-1)
- Sort blocks: if (p,q) > (r,s) lexicographically, swap blocks (no sign change)
- This normalizes all equivalent ERI representations to a canonical form
- After canonicalization, tensor list is sorted (products are commutative)
- Terms with identical canonical forms merge and sum their factors

---

## Testing Strategy

### Why gradual "swap in" is limited

pdaggerq's Python API exposes only `strings()` — the fully processed final
output. The C++ `get_ordered_strings()` method (which gives access to internal
`pq_string` objects with deltas, operators, etc.) is **not bound in Python**.
This means we cannot query pdaggerq's intermediate state.

Consequence: most pipeline steps (delta elimination, filter, canonicalize,
expansion, normal ordering) require intermediate forms that only exist within
our own pipeline. There is no pdaggerq oracle we can compare against mid-stream.

The only integration-level swap that was possible was `format_strings` (the
very last step), which operates on the same string data that pdaggerq emits.
That swap is already done and confirmed.

### Actual strategy

```
For each step:
  1. Write isolated test with hand-crafted inputs
  2. Run it
  3. Fix until it passes
  4. STOP. Report to user. Wait for feedback before proceeding.

Once ALL steps have passing isolated tests:
  5. Wire the full pipeline in WickHelper
  6. Replace ours_mp2() with our full WickHelper call
  7. Run test_mp2.py → must match pdaggerq strings() exactly
  8. STOP. Report to user. Wait for feedback.
  9. Repeat for test_mp3.py
```

We never advance to the next step's isolated test without user confirmation.
We never attempt the full integration swap until ALL isolated tests pass.

---

## Progress

| Step | Module | Isolated test | Status |
|------|--------|--------------|--------|
| Output formatting | `wick/output.py` | `tests/test_step2_output.py` | ✅ **DONE** (3/3 tests pass) |
| Delta elimination | `wick/deltas.py` | `tests/test_step3_deltas.py` | ✅ **DONE** (8/8 tests pass) |
| Label canonicalization | `wick/labels.py` | `tests/test_step4_labels.py` | ✅ **DONE** (6/6 tests pass) |
| Contraction filter | `wick/filter.py` | `tests/test_step5_filter.py` | ✅ **DONE** (5/5 tests pass) |
| General label expansion | `wick/expand.py` | `tests/test_step6_expand.py` | ✅ **DONE** (8/8 tests pass) |
| Cleanup / cancel terms | `wick/cleanup.py` | `tests/test_step7_cleanup.py` | ✅ **DONE** (7/7 tests pass) — ERI canonicalization added |
| Normal ordering | `wick/normal_order.py` | `tests/test_step8_normal_order.py` | ✅ **DONE** (15/15 tests pass) |
| Full integration MP2 | `wick/helper.py` | `tests/test_mp2.py` | ✅ **DONE** — output matches pdaggerq exactly |
| Full integration MP3 | `wick/helper.py` | `tests/test_mp3.py` | ✅ **DONE** — output matches pdaggerq exactly |

**Test Suite**: 54/54 tests passing ✅

---

## Current State of Integration Tests

### `ours_mp2()` and `ours_mp3()`

```python
def ours_mp2():
    # Full WickHelper pipeline now active
    wh = WickHelper()
    wh.add_term(1/16, ['g(p,q,r,s)', 'g(t,u,v,w)'], ['+p','+q','-s','-r','+t','+u','-w','-v'])
    wh.simplify()
    return wh.strings()
```

Output matches pdaggerq exactly for both MP2 and MP3.

### What makes it work

1. **expand_general_labels**: Expands `p,q,r,s` into 16 terms (all occ, all virt, mixed)
2. **normal_order_fermi_vacuum**: Applies Wick's theorem, generates swap and contraction terms
3. **filter_fully_contracted**: Keeps only terms with all operators contracted
4. **apply_deltas**: Eliminates delta functions by label substitution
5. **canonicalize_labels**: Renames o1,o2,... → i,j,... and v1,v2,... → a,b,...
6. **cancel_terms with ERI canonicalization**:
   - Normalizes each `<p,q||r,s>` to canonical form
   - Sorts tensors in each term (commutative products)
   - Groups by canonical form and sums factors
7. **format_strings**: Converts to pdaggerq output format

---

## Key Design Decisions (do not relitigate)

- Pure Python, Fermi vacuum only
- Antisymmetrized two-electron integrals `<p,q||r,s>`
- No one-body `f` operator for MP2/MP3 (HF reference, f contributions vanish)
- No BCH/commutator input API (we take raw operator strings)
- Output format must match pdaggerq `strings()` exactly for string comparison
- No spin blocking, no pq_graph, no code generation — out of scope
- ERI canonicalization: lexicographic sorting of bra and ket pairs, with sign tracking
  - This is essential to group equivalent terms under permutations

## Completed Milestones

✅ **Full pipeline implementation complete**
- All 7 processing steps working correctly
- All 54 unit tests passing
- MP2 and MP3 integration tests matching pdaggerq exactly
- Ready for production use or further extensions
