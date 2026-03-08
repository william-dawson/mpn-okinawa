//! Step 4b: Reclassify fluctuation-potential one-body pieces.
//!
//! Convert occ_repulsion tensors into equivalent ERI tensors:
//!   <?i||?j>  ->  <i,k||j,k>

use std::collections::HashSet;

use crate::term::{Tensor, Term};

const OCC_PREF: &[&str] = &["i", "j", "k", "l", "m", "n"];

fn all_labels(term: &Term) -> HashSet<String> {
    let mut labels = HashSet::new();
    for t in &term.tensors {
        for idx in &t.indices {
            labels.insert(idx.clone());
        }
    }
    for d in &term.deltas {
        labels.insert(d.i.clone());
        labels.insert(d.j.clone());
    }
    for op in &term.operators {
        labels.insert(op.label.clone());
    }
    labels
}

fn fresh_occ_label(used: &HashSet<String>) -> String {
    for lbl in OCC_PREF {
        if !used.contains(*lbl) {
            return lbl.to_string();
        }
    }
    let mut n = 0;
    loop {
        let lbl = format!("o{}", n);
        if !used.contains(&lbl) {
            return lbl;
        }
        n += 1;
    }
}

pub fn reclassify_occ_repulsion(terms: &[Term]) -> Vec<Term> {
    let mut out = Vec::new();
    for term in terms {
        let mut used = all_labels(term);
        let mut new_tensors = Vec::new();
        for t in &term.tensors {
            if t.name != "occ_repulsion" {
                new_tensors.push(t.clone());
                continue;
            }
            let i = &t.indices[0];
            let j = &t.indices[1];
            let k = fresh_occ_label(&used);
            used.insert(k.clone());
            new_tensors.push(Tensor::from_strings(
                "g".to_string(),
                vec![i.clone(), k.clone(), j.clone(), k],
            ));
        }
        out.push(Term::new(
            term.factor,
            term.sign,
            new_tensors,
            term.operators.clone(),
            term.deltas.clone(),
        ));
    }
    out
}
