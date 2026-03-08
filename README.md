# Wick's Theorem Engine

A Rust Wick-contraction pipeline for MP-style vacuum expectation values, with CLI comparison against `pdaggerq`.

## Project Status Warning

This codebase was generated and iterated during a hackathon (CompPhysHack 2026): https://qc-hybrid.github.io/CompPhysHack2026/

It is a rapid prototype and should not be treated as production-ready software.

## Building

```bash
cargo build --release
```

The binary is installed as `wick` (see `target/release/wick`).

## Attribution and Provenance

This repository was developed with LLM assistance while giving the model direct access to the `pdaggerq` source code (https://github.com/edeprince3/pdaggerq). As a result, parts of this implementation were likely heavily influenced by `pdaggerq`, and were probably plagiarized heavily from its code structure and logic.

The project is licensed under Apache License 2.0 to match `pdaggerq` (see `LICENSE`).

## Simplification Strategy

The project uses **manual consolidation rules** (`step_6_cleanup::cancel_terms`) for MP-style expressions.

## Pipeline (Default)

`WickHelper::simplify()` runs:

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

- `src/term.rs`: core data structures and index-space helpers
- `src/helper.rs`: pipeline orchestrator
- `src/step_1_expand.rs`: general-label expansion
- `src/step_2_normal_order.rs`: Wick normal-order tree generation
- `src/step_3_filter.rs`: fully-contracted filtering
- `src/step_4_deltas.rs`: delta elimination
- `src/step_4b_reclassify.rs`: occ-repulsion reclassification
- `src/step_5_labels.rs`: pdaggerq-style label canonicalization
- `src/step_6_cleanup.rs`: default manual cleanup/consolidation
- `src/step_6a_linked_cluster.rs`: linked-diagram filter
- `src/step_7_output.rs`: formatting to pdaggerq-style strings
- `src/step_8_denominator.rs`: denominator extraction
- `src/main.rs`: command-line interface

## Tests

```bash
cargo test
```
