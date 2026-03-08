# Wick's Theorem Engine

A pure Python implementation of a Wick-contraction pipeline for MP-style vacuum expectation values, with side-by-side CLI comparison against `pdaggerq`.

## Current Status

- MP1: matches `pdaggerq` in coefficient/structure up to dummy-index orientation.
- MP2: parity-mode output has the correct two-term structure (orientation differences remain).
- MP3: currently gives 5 terms in parity mode; 3 match textbook topologies (pp/hh/ph), 2 are disconnected leftovers.

Important: this repository has two useful output modes for MP examples:

- **Reference-parity mode** (`filter_unlinked=False`): keeps disconnected terms, closer to raw `pdaggerq strings()` behavior.
- **Linked-energy mode** (`filter_unlinked=True`, default): applies linked-diagram filtering to show the connected contribution typically used in final MP energy expressions.

## CLI Comparison (pdaggerq first, ours second)

The CLI now prints both outputs so you can inspect differences directly:

```bash
python wick_cli.py --mp2
```

Order of display:
1. `pdaggerq` reference `strings()` output.
2. Our output from `WickHelper`.
3. `Exact string match vs pdaggerq: ...`

## Pipeline (As Implemented)

`WickHelper.simplify()` currently runs:

1. `expand_general_labels`
2. `normal_order_fermi_vacuum`
3. `filter_fully_contracted`
4. `apply_deltas`
5. `reclassify_occ_repulsion`
6. `canonicalize_labels_pdaggerq_style`
7. `cancel_terms`
8. `filter_energy_subspace` (optional)
9. `filter_unlinked_diagrams` (optional)

## MP2 Walkthrough: Exact Current Behavior

Input used in CLI/tests:

- factor: `1/16`
- tensors: `['g(p,q,r,s)', 'g(t,u,v,w)']`
- ops: `['+p','+q','-r','-s','+t','+u','-v','-w']`

### Step 0: Input

- Terms: `1`

### Step 1: General Label Expansion

- Function: `expand_general_labels`
- Output terms: `256`
- Reason: 8 general labels -> `2^8` occupied/virtual assignments.

### Step 2: Fermi-Vacuum Normal Ordering

- Function: `normal_order_fermi_vacuum`
- Output terms: `2840`
- Behavior: recursive swap/contraction branching.

### Step 3: Keep Fully Contracted Terms

- Function: `filter_fully_contracted`
- Output terms: `24`

### Step 4: Delta Elimination

- Function: `apply_deltas`
- Output terms: `24`

### Step 5: Reclassify occ_repulsion

- Function: `reclassify_occ_repulsion`
- Output terms: `24`

### Step 6: Canonical Label Rename

- Function: `canonicalize_labels_pdaggerq_style`
- Output terms: `24`

### Step 7: Cleanup / Merge Equivalent Terms

- Function: `cancel_terms`
- Output terms: `3`
- Terms:
  - `['+0.250', '<a,b||i,j>', '<i,j||a,b>']`
  - `['+0.250', '<i,j||i,j>', '<k,l||k,l>']`
  - `['-1.00', '<a,i||i,j>', '<j,k||a,k>']`

### Step 8: Energy-Subspace Projection

- Function: `filter_energy_subspace`
- Output terms: `2`
- Removed: odd occupied/virtual balance term (singles-like residual)
- Remaining:
  - `['+0.250', '<a,b||i,j>', '<i,j||a,b>']`
  - `['+0.250', '<i,j||i,j>', '<k,l||k,l>']`

### Step 9: Linked-Diagram Filter (default ON)

- Function: `filter_unlinked_diagrams`
- Output terms: `1`
- Remaining linked term:
  - `['+0.250', '<a,b||i,j>', '<i,j||a,b>']`

This is the default `WickHelper()` MP2 result.

## Why MP2 Can Show One or Two Terms

For MP2 in this codebase:

- `filter_unlinked=False` (parity-style view): two terms survive after projection.
- `filter_unlinked=True` (default linked-energy view): only the connected doubles term remains.

So seeing one vs two terms depends on whether you are looking at disconnected/reference pieces or only linked contributions.

## MP3: What We Have Right Now

In parity mode (`filter_unlinked=False`, `project_energy_subspace=False`), current output has 5 terms:

- `['+0.1250', '<a,b||c,d>', '<c,d||i,j>', '<i,j||a,b>']`  -> pp topology
- `['+0.1250', '<a,b||i,j>', '<i,j||k,l>', '<k,l||a,b>']`  -> hh topology
- `['-1.00', '<a,c||i,k>', '<b,k||c,j>', '<i,j||a,b>']`    -> ph topology
- `['-0.1250', '<i,j||i,j>', '<k,l||k,l>', '<m,n||m,n>']`  -> disconnected triple bubble
- `['-0.3750', '<a,b||i,j>', '<i,j||a,b>', '<k,l||k,l>']`  -> semi-disconnected bubble product

Interpretation:

- The first 3 correspond to textbook MP3 classes (pp/hh/ph).
- The last 2 are disconnected leftovers that should cancel/drop for connected MP3 correlation energy.

## Modules

- `wick/term.py`: core data structures and index-space helpers
- `wick/helper.py`: pipeline orchestrator
- `wick/step_1_expand.py`: general-label expansion
- `wick/step_2_normal_order.py`: Wick normal-order tree generation
- `wick/step_3_filter.py`: fully-contracted filtering
- `wick/step_4_deltas.py`: delta elimination
- `wick/step_4b_reclassify.py`: occ-repulsion reclassification
- `wick/step_5_labels.py`: pdaggerq-style label canonicalization
- `wick/step_6_cleanup.py`: integral canonicalization + term merging
- `wick/step_6b_projection.py`: energy-subspace projection
- `wick/step_6a_linked_cluster.py`: linked-diagram filter
- `wick/step_7_output.py`: formatting to pdaggerq-style strings
- `wick/step_8_denominator.py`: denominator extraction
- `wick_cli.py`: interactive/non-interactive CLI with reference comparison

## Tests

Run targeted MP checks:

```bash
pytest -q tests/test_mp1.py tests/test_mp2.py tests/test_mp3.py
```

Run full local tests:

```bash
pytest -q tests/
```

## Open Work

Primary remaining task is to eliminate the two disconnected MP3 leftovers by fixing earlier cancellation/consolidation so connected MP3 reduces cleanly to the expected pp/hh/ph content.
