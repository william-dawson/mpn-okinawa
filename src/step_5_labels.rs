//! Step 5: Label Canonicalization
//!
//! Rename summed index labels to conventional names:
//!   occupied: i, j, k, l, m, n
//!   virtual:  a, b, c, d, e, f

use std::collections::{HashMap, HashSet};

use crate::term::Term;

fn collect_labels_in_order(term: &Term) -> Vec<String> {
    let mut seen = HashSet::new();
    let mut order = Vec::new();

    let mut visit = |lbl: &str| {
        if seen.insert(lbl.to_string()) {
            order.push(lbl.to_string());
        }
    };

    for t in &term.tensors {
        for idx in &t.indices {
            visit(idx);
        }
    }
    for d in &term.deltas {
        visit(&d.i);
        visit(&d.j);
    }
    for op in &term.operators {
        visit(&op.label);
    }

    order
}

pub fn canonicalize_labels_pdaggerq_style(terms: &[Term]) -> Vec<Term> {
    let occ_in: Vec<String> = (0..30).map(|i| format!("o{}", i)).collect();
    let virt_in: Vec<String> = (0..30).map(|i| format!("v{}", i)).collect();
    let occ_out: Vec<String> = "ijklmn".chars().map(|c| c.to_string()).collect();
    let virt_out: Vec<String> = "abcdef".chars().map(|c| c.to_string()).collect();

    let mut out = Vec::new();
    for term in terms {
        let labels = collect_labels_in_order(term);
        let mut used: HashSet<String> = labels.iter().cloned().collect();
        let mut mapping = HashMap::new();
        let mut mapped_values: HashSet<String> = HashSet::new();

        for src in &occ_in {
            if !used.contains(src) {
                continue;
            }
            for dst in &occ_out {
                if !used.contains(dst) && !mapped_values.contains(dst) {
                    mapping.insert(src.clone(), dst.clone());
                    mapped_values.insert(dst.clone());
                    used.insert(dst.clone());
                    break;
                }
            }
        }

        for src in &virt_in {
            if !used.contains(src) {
                continue;
            }
            for dst in &virt_out {
                if !used.contains(dst) && !mapped_values.contains(dst) {
                    mapping.insert(src.clone(), dst.clone());
                    mapped_values.insert(dst.clone());
                    used.insert(dst.clone());
                    break;
                }
            }
        }

        out.push(term.renamed(&mapping));
    }
    out
}
