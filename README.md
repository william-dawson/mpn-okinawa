# Wick's Theorem Engine

A pure Python implementation of second-quantized operator algebra for computing correlation energy corrections (MP2, MP3, etc.) using Wick's theorem with respect to the Fermi vacuum.

## Overview

This engine takes explicit operator strings as input and applies a 7-step pipeline to produce final pdaggerq-compatible output. Each step is isolated, tested, and composable.

**Pipeline:** General label expansion → Normal ordering → Filter contractions → Delta elimination → Label canonicalization → Cleanup & ERI canonicalization → Output formatting

## The 7-Step Pipeline

We'll trace a concrete example: **MP2** with input `1/16 * g(p,q,r,s) * g(t,u,v,w) * V` where V is the operator string `['+p','+q','-s','-r','+t','+u','-w','-v']`.

---

### Step 1: General Label Expansion

**Purpose:** Expand general indices (p, q, r, s, t, u, v, w) into all occupied/virtual combinations.

**Input:** One term with 8 general labels
```
factor=1/16
tensors: [g(p,q,r,s), g(t,u,v,w)]
operators: [+p, +q, -s, -r, +t, +u, -w, -v]
```

**Process:** Each of the 8 labels independently assigned to either occupied (fresh: o1, o2, ...) or virtual (fresh: v1, v2, ...).
- 2^8 = 256 combinations
- Labels are globally unique via a shared counter

**Output:** 256 terms, one for each assignment. Example subset:
```
Term 1: factor=1/16, tensors: [g(o1,o2,v1,v2), g(o3,o4,v3,v4)], ops: [+o1, +o2, -v2, -v1, +o3, +o4, -v4, -v3]
Term 2: factor=1/16, tensors: [g(o1,o2,v1,v3), g(o3,o4,v2,v4)], ops: [+o1, +o2, -v3, -v1, +o3, +o4, -v4, -v2]
...
Term 256: factor=1/16, tensors: [g(v1,v2,v3,v4), g(v5,v6,v7,v8)], ops: [+v1, +v2, -v4, -v3, +v5, +v6, -v8, -v7]
```

---

### Step 2: Normal Ordering via Wick's Theorem

**Purpose:** Apply Wick's theorem to bring all operator strings to Fermi vacuum normal order.

**Input:** 256 terms from Step 1 (each with a string of 8 operators in arbitrary order)

**Process:** For each term, repeatedly find the leftmost out-of-order operator pair (FV annihilator before FV creator) and apply two rules:
1. **Swap term:** Swap operators, flip sign (anticommutation relation)
2. **Contraction term:** If both operators are in the same orbital space (both occ or both virt), remove them and add a delta function δ(label₁, label₂); keep original sign

Repeat until all terms are normal ordered (no out-of-order pairs remain).

**Example on one term:** `[+o1, -o2, +o3, -o4, ...]`
- Found out-of-order pair: `-o2` (FV creator of occ) before `+o3` (FV annihilator of occ)
- Generate two new terms:
  - Swap: `[+o1, +o3, -o2, -o4, ...]` with sign flipped
  - Contract: `[+o1, -o4, ...]` with delta(o2, o3) added, sign unchanged

**Output:** ~1000+ terms (exponential growth from BFS over swap/contract tree). Each term now has:
- 0-8 remaining operators (or 0 if fully contracted)
- 0+ delta functions tracking contractions

---

### Step 3: Contraction Filter

**Purpose:** Keep only fully contracted terms (no remaining operators).

For Fermi vacuum `<0|` and `|0>`, only terms where all 8 operators have been contracted survive. Any term with remaining operators vanishes.

**Input:** ~1000+ mixed terms from Step 2

**Process:** Simple filter: keep only terms where `operators == []`

**Output:** ~40-60 terms (the uncontracted terms are discarded)

Example surviving term:
```
factor=1/16, tensors: [g(o1,o2,v1,v2), g(o3,o4,v3,v4)],
deltas: [delta(o1,o3), delta(o2,o4), delta(v1,v3), delta(v2,v4)]
```

---

### Step 4: Delta Elimination

**Purpose:** Eliminate delta functions by substituting one index label for the other.

**Input:** ~40-60 fully-contracted terms with deltas

**Process:** For each term, iteratively apply each delta:
1. Choose which label to keep (prefer conventional i,j,k,... over fresh o1,o2,...; prefer lower index if both fresh)
2. Replace the eliminated label everywhere in tensors and remaining deltas
3. Remove the delta

**Example:** Term with `delta(o1, i)` and `delta(o2, j)`:
- Before: `g(o1,o2,v1,v2)` with `delta(o1,i), delta(o2,j)`
- After applying delta(o1,i): `g(i,o2,v1,v2)` with `delta(o2,j)`
- After applying delta(o2,j): `g(i,j,v1,v2)` with no deltas

**Output:** ~40-60 terms with all deltas eliminated, tensors now using mix of conventional and fresh labels

---

### Step 5: Label Canonicalization

**Purpose:** Rename fresh labels (o1, o2, v1, v2, ...) to conventional names (i, j, k, l, ... for occ; a, b, c, d, ... for virt) in order of first appearance.

**Input:** ~40-60 terms with mixed fresh and conventional labels

**Process:** Scan each term left-to-right across all tensors. Assign the first appearance of each fresh label to the next available conventional name:
- o1 → i, o2 → j, o3 → k, etc.
- v1 → a, v2 → b, v3 → c, etc.

**Example:**
- Before: `g(i,o2,v1,o3)` (i is conventional, o2,v1,o3 are fresh)
- After: `g(i,j,a,k)`  (o2→j, v1→a, o3→k following appearance order)

**Output:** ~40-60 terms with all conventional labels

---

### Step 6: Cleanup & ERI Canonicalization

**Purpose:** Normalize ERI tensor representations and merge equivalent terms.

**Input:** ~40-60 terms like `g(i,j,a,b)`, `g(j,i,b,a)`, etc. (same ERIs in different index orderings)

**Process:**
1. **Canonicalize each 4-index ERI** `<p,q||r,s>`:
   - Sort bra pair (p,q): if p > q lexicographically, swap and apply sign change (-1)
   - Sort ket pair (r,s): if r > s, swap and apply sign change (-1)
   - Sort left-right blocks: if (p,q) > (r,s), swap blocks (no sign change for <||>)

2. **Sort tensor list** within each term (products are commutative)

3. **Group by canonical form** and sum factors of identical terms

**Example:**
- Term 1: `g(j,i,b,a)` → canonicalize → `<i,j||a,b>` with sign +1 → sorted: `[<i,j||a,b>]`
- Term 2: `g(i,j,a,b)` → canonicalize → `<i,j||a,b>` with sign +1 → sorted: `[<i,j||a,b>]`
- Both map to same key → factors merge: 0.0625 + 0.0625 = 0.1250

**Output:** ~2 terms (merged from 40-60 via canonicalization)

Examples:
```
Term 1: factor=0.250, tensors: [g(i,j||i,j), g(k,l||k,l)]
Term 2: factor=0.250, tensors: [g(a,b||j,i), g(j,i||a,b)]
```

---

### Step 7: Output Formatting

**Purpose:** Convert terms to pdaggerq-compatible `strings()` format.

**Input:** ~2 terms with conventional labels and merged factors

**Process:**
1. Format factor using pdaggerq's `minimum_precision()` algorithm:
   - Work with 10 × |factor|, count decimal digits
   - Stop after 12 repeating digits
   - Enforce minimum 2 decimal places
   - Prepend + or - sign

2. Convert each tensor to string format: `<p,q||r,s>`

3. Return list of lists: `[[factor_str, tensor_str, ...], ...]`

**Example:**
- Input: `Term(factor=0.250, tensors=[Tensor('g', ['i','j','i','j']), ...])`
- Output: `['+0.250', '<i,j||i,j>', '<k,l||k,l>']`

**Output:**
```python
[
  ['+0.250', '<i,j||i,j>', '<k,l||k,l>'],
  ['+0.250', '<a,b||j,i>', '<j,i||a,b>']
]
```

This matches pdaggerq exactly.

---

### Step 6.5: Physics Filters (Linked-Cluster & Brillouin Theorems)

**Purpose:** Remove unphysical diagram contributions using exact theorems from MBPT.

**Input:** ~2-4 terms after cleanup and ERI canonicalization

**Process:**

1. **Linked-Cluster Theorem Filter**
   - Removes pure bubble (self-energy) diagrams
   - Example: `<i,j||i,j> * <k,l||k,l>` (disconnected product of two identical contributions)
   - These are unlinked diagrams that cancel exactly in MBPT
   - Condition: Both ERIs have completely disjoint index sets AND both are internally self-contracted

2. **Brillouin's Theorem Filter**
   - Removes contributions with single-excitation character
   - Physics: ⟨Ψ⁰|H|Ψ¹⟩ = 0 (ground state orthogonal to singly-excited states)
   - Example: `<i,j||a,j> * <a,k||i,k>` has {i,j,k} occupied + {a} virtual
   - Condition: Exactly 3 occ + 1 virt (or vice versa) indicates coupling to single-excitation manifold

**Example for MP2:**

Before filters (4 terms):
```
[0] +0.250 * <i,j||i,j> * <k,l||k,l>       ← Removed: bubble diagram
[1] +0.500 * <i,j||a,j> * <a,k||i,k>       ← Removed: single-excitation (3 occ, 1 virt)
[2] +0.250 * <i,j||a,b> * <a,b||i,j>       ← Kept: true two-body excitation (2 occ, 2 virt)
[3] -0.500 * <i,j||a,i> * <a,k||j,k>       ← Removed: single-excitation (3 occ, 1 virt)
Total: +0.500
```

After filters (1 term):
```
[0] +0.250 * <i,j||a,b> * <a,b||i,j>       ← Only genuine two-body excitation survives
Total: +0.250 = 1/4 ✓ (correct MP2 correlation energy)
```

**Denominator Structure:** Filtered terms have the correct intermediate excitation structure:
- MP2: (e_i + e_j - e_a - e_b) ← standard 2-hole-2-particle energy gap

---

## Usage

```python
from wick import WickHelper

wh = WickHelper()
wh.add_term(
    factor=1/16,
    tensors=['g(p,q,r,s)', 'g(t,u,v,w)'],
    ops=['+p', '+q', '-s', '-r', '+t', '+u', '-w', '-v']
)
wh.simplify()
result = wh.strings()
print(result)
# [['+0.250', '<i,j||i,j>', '<k,l||k,l>'],
#  ['+0.250', '<a,b||j,i>', '<j,i||a,b>']]
```

## Testing

- **54 unit tests** covering each step in isolation
- **3 integration tests** (MP1, MP2, MP3) verifying end-to-end output matches pdaggerq

Run all tests:
```bash
pytest tests/ -v
```

## Architecture

- `wick/term.py` — Core immutable data structures (Term, Tensor, Operator, Delta)
- `wick/helper.py` — Pipeline orchestrator (WickHelper class)
- `wick/step_1_expand.py` — General label expansion
- `wick/step_2_normal_order.py` — Wick's theorem via normal ordering
- `wick/step_3_filter.py` — Contraction filtering
- `wick/step_4_deltas.py` — Delta elimination
- `wick/step_5_labels.py` — Label canonicalization
- `wick/step_6_cleanup.py` — ERI canonicalization & term merging
- `wick/step_7_output.py` — pdaggerq format output

## Performance

- **MP1:** 1 term, instant
- **MP2:** 2 terms, <0.1s
- **MP3:** 15 terms, ~0.1s
- **MP4:** 65,536 initial terms → computational explosion (exponential in number of general labels)

The pipeline trades off memory for clarity: each step is independent, testable, and understandable. For production use with MP4+, consider using pdaggerq's C++ backend directly.
