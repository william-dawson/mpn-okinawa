# Wick's Theorem Engine

A pure Python Wick-contraction pipeline for MP-style vacuum expectation values, with CLI comparison against `pdaggerq`.

## Installation

Base package:

```bash
pip install -e .
```

Test/development extras (includes optional `pdaggerq` from GitHub):

```bash
pip install -e '.[test]'
```

## Simplification Strategy

The project uses **manual consolidation rules** (`step_6_cleanup.cancel_terms`) for MP-style expressions.

## Pipeline (Default)

`WickHelper.simplify()` runs:

1. `expand_general_labels`
2. `normal_order_fermi_vacuum`
3. `filter_fully_contracted`
4. `apply_deltas`
5. `reclassify_occ_repulsion`
6. `canonicalize_labels_pdaggerq_style`
7. `cancel_terms` (default manual consolidation)
8. `filter_unlinked_diagrams` (optional)

## MP2 Workflow (Theory + Verbose Trace)

You can regenerate this walkthrough with:

```bash
wick --mp2 --verbose
```

For MP2 (`['v','v']`), the expected role of each step and the current observed output are:

| Step | Stage | Theory | Observed (`wick --mp2 --verbose`) |
|---|---|---|---|
| 1 | `expand_general_labels` | Expand each general index into occupied/virtual assignments. | `400` terms |
| 2 | `normal_order_fermi_vacuum` | Apply Wick reordering and contraction branching. | `3446` terms |
| 3 | `filter_fully_contracted` | Keep only vacuum-surviving fully contracted terms. | `38` terms |
| 4 | `apply_deltas` | Eliminate Kronecker deltas by substitution. | `38` terms (structure simplified) |
| 5 | `reclassify_occ_repulsion` | Convert `occ_repulsion` intermediates into ERIs. | `38` terms (ERI form) |
| 6 | `canonicalize_labels_pdaggerq_style` | Normalize dummy labels to conventional symbols. | `38` terms |
| 7 | `cancel_terms` | Merge equivalent tensor products and sum coefficients. | `2` terms |
| 8 | `filter_unlinked_diagrams` | Remove disconnected contributions; keep linked energy terms. | `1` term |

Final MP2 term:

```text
+0.250 * <a,b||i,j> * <i,j||a,b> / (e_i + e_j - e_a - e_b)
```

## CLI

Run MP examples:

```bash
wick --mp1
wick --mp2
wick --mp3
```

## Modules

- `wick/term.py`: core data structures and index-space helpers
- `wick/helper.py`: pipeline orchestrator
- `wick/step_1_expand.py`: general-label expansion
- `wick/step_2_normal_order.py`: Wick normal-order tree generation
- `wick/step_3_filter.py`: fully-contracted filtering
- `wick/step_4_deltas.py`: delta elimination
- `wick/step_4b_reclassify.py`: occ-repulsion reclassification
- `wick/step_5_labels.py`: pdaggerq-style label canonicalization
- `wick/step_6_cleanup.py`: default manual cleanup/consolidation
- `wick/step_6a_linked_cluster.py`: linked-diagram filter
- `wick/step_7_output.py`: formatting to pdaggerq-style strings
- `wick/step_8_denominator.py`: denominator extraction
- `wick/cli.py`: interactive/non-interactive CLI

## Tests

```bash
pytest -q tests
```
