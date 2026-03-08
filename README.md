# Wick's Theorem Engine

A pure Python Wick-contraction pipeline for MP-style vacuum expectation values, with CLI comparison against `pdaggerq`.

## Simplification Strategy

The project currently uses **manual consolidation rules by default** (`step_6_cleanup.cancel_terms`) because this path is currently the best balance of correctness and speed for MP1/MP2/MP3 runs.

A graph-isomorphism consolidator is available as an **experimental alternative** (`step_6d_graph_isomorphism.cancel_terms_graph_isomorphic`) and can be enabled from the CLI.

## Pipeline (Default)

`WickHelper.simplify()` runs:

1. `expand_general_labels`
2. `normal_order_fermi_vacuum`
3. `filter_fully_contracted`
4. `apply_deltas`
5. `reclassify_occ_repulsion`
6. `canonicalize_labels_pdaggerq_style`
7. `cancel_terms` (default manual consolidation)
8. `filter_connected_terms_graph` (optional)
9. `filter_unlinked_diagrams` (optional)

## CLI

Run MP examples:

```bash
python wick_cli.py --mp1
python wick_cli.py --mp2
python wick_cli.py --mp3
```

Try graph-isomorphic simplification (experimental):

```bash
python wick_cli.py --mp3 --graph-simplify
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
- `wick/step_6d_graph_isomorphism.py`: experimental graph-isomorphic consolidation
- `wick/step_6c_graph_connected.py`: optional connected-term filter
- `wick/step_6a_linked_cluster.py`: linked-diagram filter
- `wick/step_7_output.py`: formatting to pdaggerq-style strings
- `wick/step_8_denominator.py`: denominator extraction
- `wick_cli.py`: interactive/non-interactive CLI

## Tests

```bash
pytest -q tests
```
