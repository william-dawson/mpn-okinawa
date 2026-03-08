//! Step 2: Normal Ordering via Wick's Theorem (Fermi Vacuum)
//!
//! Repeatedly scan for the leftmost pair of adjacent operators that is
//! out of Fermi-vacuum normal order (FV annihilator before FV creator).
//! Generate swap + contraction terms until all terms are normal ordered.

use std::collections::VecDeque;

use crate::term::{is_occ, is_virt, Delta, Operator, Term};

fn first_out_of_order(operators: &[Operator]) -> Option<usize> {
    for i in 0..operators.len().saturating_sub(1) {
        if !operators[i].dagger_fermi() && operators[i + 1].dagger_fermi() {
            return Some(i);
        }
    }
    None
}

fn swap_one(term: &Term) -> Vec<Term> {
    let ops = &term.operators;
    let i = match first_out_of_order(ops) {
        Some(i) => i,
        None => return vec![term.clone()],
    };

    let ann = &ops[i];
    let cre = &ops[i + 1];

    // 1. Swap term
    let mut swap_ops: Vec<Operator> = ops[..i].to_vec();
    swap_ops.push(cre.clone());
    swap_ops.push(ann.clone());
    swap_ops.extend_from_slice(&ops[i + 2..]);

    let swap_term = Term::new(
        term.factor,
        -term.sign,
        term.tensors.clone(),
        swap_ops,
        term.deltas.clone(),
    );

    let mut results = vec![swap_term];

    // 2. Contraction term (same space only)
    let both_occ = is_occ(&ann.label) && is_occ(&cre.label);
    let both_virt = is_virt(&ann.label) && is_virt(&cre.label);

    if both_occ || both_virt {
        let mut contr_ops: Vec<Operator> = ops[..i].to_vec();
        contr_ops.extend_from_slice(&ops[i + 2..]);

        let (d1, d2) = if ann.label <= cre.label {
            (ann.label.clone(), cre.label.clone())
        } else {
            (cre.label.clone(), ann.label.clone())
        };

        let mut contr_deltas = term.deltas.clone();
        contr_deltas.push(Delta { i: d1, j: d2 });

        let contr_term = Term::new(
            term.factor,
            term.sign,
            term.tensors.clone(),
            contr_ops,
            contr_deltas,
        );
        results.push(contr_term);
    }

    results
}

pub fn normal_order_fermi_vacuum(terms: &[Term]) -> Vec<Term> {
    let mut ordered = Vec::new();
    let mut queue: VecDeque<Term> = terms.iter().cloned().collect();

    while let Some(term) = queue.pop_front() {
        if term.is_normal_ordered_fermi() {
            ordered.push(term);
        } else {
            for t in swap_one(&term) {
                queue.push_back(t);
            }
        }
    }

    ordered
}
