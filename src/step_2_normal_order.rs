//! Step 2: Normal Ordering via Wick's Theorem (Fermi Vacuum)
//!
//! Repeatedly scan for the leftmost pair of adjacent operators that is
//! out of Fermi-vacuum normal order (FV annihilator before FV creator).
//! Generate swap + contraction terms until all terms are normal ordered.

use rustc_hash::FxHashMap;

use crate::term::{is_occ, is_virt, Delta, Operator, Term};

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum Space {
    Occ,
    Virt,
    Other,
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct IntOperator {
    label: u16,
    dagger: bool,
    space: Space,
}

#[derive(Clone, Copy)]
struct IntDelta {
    i: u16,
    j: u16,
}

fn dagger_fermi(op: &IntOperator) -> bool {
    match op.space {
        Space::Occ => !op.dagger,
        Space::Virt => op.dagger,
        Space::Other => panic!(
            "dagger_fermi undefined for general labels in normal ordering; expand labels first"
        ),
    }
}

fn first_out_of_order(operators: &[IntOperator]) -> Option<usize> {
    for i in 0..operators.len().saturating_sub(1) {
        if !dagger_fermi(&operators[i]) && dagger_fermi(&operators[i + 1]) {
            return Some(i);
        }
    }
    None
}

#[derive(Clone)]
struct Expansion {
    rel_sign: i32,
    rel_deltas: Vec<IntDelta>,
    operators: Vec<IntOperator>,
}

fn expand_ops_memo(
    ops: &[IntOperator],
    cache: &mut FxHashMap<Vec<IntOperator>, usize>,
    memo_results: &mut Vec<Vec<Expansion>>,
) -> usize {
    if let Some(&hit) = cache.get(ops) {
        return hit;
    }

    let result = match first_out_of_order(ops) {
        None => vec![Expansion {
            rel_sign: 1,
            rel_deltas: Vec::new(),
            operators: ops.to_vec(),
        }],
        Some(i) => {
            let ann = &ops[i];
            let cre = &ops[i + 1];

            let mut out = Vec::new();

            // 1. Swap branch: sign flip.
            let mut swap_ops: Vec<IntOperator> = ops[..i].to_vec();
            swap_ops.push(*cre);
            swap_ops.push(*ann);
            swap_ops.extend_from_slice(&ops[i + 2..]);
            let swap_children_idx = expand_ops_memo(&swap_ops, cache, memo_results);
            for child in memo_results[swap_children_idx].iter() {
                out.push(Expansion {
                    rel_sign: -child.rel_sign,
                    rel_deltas: child.rel_deltas.clone(),
                    operators: child.operators.clone(),
                });
            }

            // 2. Contraction branch (same space only).
            let both_occ = ann.space == Space::Occ && cre.space == Space::Occ;
            let both_virt = ann.space == Space::Virt && cre.space == Space::Virt;
            if both_occ || both_virt {
                let mut contr_ops: Vec<IntOperator> = ops[..i].to_vec();
                contr_ops.extend_from_slice(&ops[i + 2..]);

                let (d1, d2) = if ann.label <= cre.label {
                    (ann.label, cre.label)
                } else {
                    (cre.label, ann.label)
                };
                let this_delta = IntDelta { i: d1, j: d2 };

                let contr_children_idx = expand_ops_memo(&contr_ops, cache, memo_results);
                for child in memo_results[contr_children_idx].iter() {
                    // Parent contraction happens earlier than child contractions.
                    let mut deltas = Vec::with_capacity(child.rel_deltas.len() + 1);
                    deltas.push(this_delta);
                    deltas.extend(child.rel_deltas.iter().cloned());
                    out.push(Expansion {
                        rel_sign: child.rel_sign,
                        rel_deltas: deltas,
                        operators: child.operators.clone(),
                    });
                }
            }

            out
        }
    };

    let idx = memo_results.len();
    memo_results.push(result);
    cache.insert(ops.to_vec(), idx);
    idx
}

fn canonicalize_ops_for_cache(ops: &[Operator]) -> (Vec<IntOperator>, Vec<u16>, Vec<String>) {
    let mut actual_label_to_id: FxHashMap<&str, u16> = FxHashMap::default();
    let mut id_to_actual_label: Vec<String> = Vec::new();
    let mut actual_ops: Vec<IntOperator> = Vec::with_capacity(ops.len());

    for op in ops {
        let id = if let Some(id) = actual_label_to_id.get(op.label.as_str()) {
            *id
        } else {
            let id = id_to_actual_label.len() as u16;
            actual_label_to_id.insert(op.label.as_str(), id);
            id_to_actual_label.push(op.label.clone());
            id
        };
        let space = if is_occ(&op.label) {
            Space::Occ
        } else if is_virt(&op.label) {
            Space::Virt
        } else {
            Space::Other
        };
        actual_ops.push(IntOperator {
            label: id,
            dagger: op.dagger,
            space,
        });
    }

    let mut actual_to_canon: Vec<u16> = vec![u16::MAX; id_to_actual_label.len()];
    let mut canon_to_actual: Vec<u16> = Vec::new();
    let mut canon_ops: Vec<IntOperator> = Vec::with_capacity(actual_ops.len());

    for op in actual_ops {
        let canon_id = if actual_to_canon[op.label as usize] != u16::MAX {
            actual_to_canon[op.label as usize]
        } else {
            let cid = canon_to_actual.len() as u16;
            actual_to_canon[op.label as usize] = cid;
            canon_to_actual.push(op.label);
            cid
        };
        canon_ops.push(IntOperator {
            label: canon_id,
            dagger: op.dagger,
            space: op.space,
        });
    }

    (canon_ops, canon_to_actual, id_to_actual_label)
}

fn remap_expansion_to_actual(
    exp: &Expansion,
    canon_to_actual: &[u16],
    id_to_actual_label: &[String],
) -> (Vec<Operator>, Vec<Delta>) {
    let mut operators = Vec::with_capacity(exp.operators.len());
    for op in exp.operators.iter() {
        operators.push(Operator {
            label: id_to_actual_label[canon_to_actual[op.label as usize] as usize].clone(),
            dagger: op.dagger,
        });
    }

    let mut deltas = Vec::with_capacity(exp.rel_deltas.len());
    for d in exp.rel_deltas.iter() {
        deltas.push(Delta {
            i: id_to_actual_label[canon_to_actual[d.i as usize] as usize].clone(),
            j: id_to_actual_label[canon_to_actual[d.j as usize] as usize].clone(),
        });
    }

    (operators, deltas)
}

pub fn normal_order_fermi_vacuum(terms: &[Term]) -> Vec<Term> {
    let mut ordered = Vec::new();
    let mut cache: FxHashMap<Vec<IntOperator>, usize> = FxHashMap::default();
    let mut memo_results: Vec<Vec<Expansion>> = Vec::new();

    for term in terms {
        let (canon_ops, canon_to_actual, id_to_actual_label) =
            canonicalize_ops_for_cache(&term.operators);
        let expansions_idx = expand_ops_memo(&canon_ops, &mut cache, &mut memo_results);
        for exp in memo_results[expansions_idx].iter() {
            let (operators, rel_deltas_actual) =
                remap_expansion_to_actual(exp, &canon_to_actual, &id_to_actual_label);
            let mut deltas = term.deltas.clone();
            deltas.extend(rel_deltas_actual);
            ordered.push(Term::new(
                term.factor,
                term.sign * exp.rel_sign,
                term.tensors.clone(),
                operators,
                deltas,
            ));
        }
    }

    ordered
}
