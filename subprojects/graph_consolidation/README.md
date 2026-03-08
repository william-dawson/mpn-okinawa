# Graph Consolidation Prototype

This subproject is an experiment toward graph/topology-based consolidation of MBPT terms.

## Goal

Provide a more physics-structured view than raw string equality by:

- Building a contraction graph for each term.
- Detecting connected vs disconnected contributions.
- Classifying MP3-like terms into `pp`, `hh`, `ph`, or `disconnected`.
- Producing a canonical graph-style key for grouping attempts.

## Files

- `prototype.py`: parser + graph builder + topology classifier + grouping helpers.
- `demo_mp3.py`: compares `pdaggerq` MP3 vs current `WickHelper` MP3 using these tools.

## Run

From project root:

```bash
python subprojects/graph_consolidation/demo_mp3.py
```

## Notes

- Canonical graph key is currently WL-style refinement, intended as a prototype signal, not a guaranteed complete isomorphism certificate.
- This is intentionally isolated from the production simplification path.
