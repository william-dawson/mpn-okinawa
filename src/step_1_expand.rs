//! Step 1: General Label Expansion
//!
//! Expand operators with general indices (p, q, r, s, ...) into
//! occupied + virtual cases. A term with n general labels expands
//! into 2^n terms.

use std::collections::HashMap;

use crate::term::{is_general, Term};

fn general_labels_in_term(term: &Term) -> Vec<String> {
    let mut seen = Vec::new();
    for op in &term.operators {
        if is_general(&op.label) && !seen.contains(&op.label) {
            seen.push(op.label.clone());
        }
    }
    for t in &term.tensors {
        for idx in &t.indices {
            if is_general(idx) && !seen.contains(idx) {
                seen.push(idx.clone());
            }
        }
    }
    seen
}

pub fn expand_general_labels(terms: &[Term], counter: &mut usize) -> Vec<Term> {
    let mut result = Vec::new();

    for term in terms {
        let gen_labels = general_labels_in_term(term);
        if gen_labels.is_empty() {
            result.push(term.clone());
            continue;
        }

        let n = gen_labels.len();
        // Iterate over all 2^n assignments
        for bits in 0..(1u32 << n) {
            let mut mapping = HashMap::new();
            for (k, label) in gen_labels.iter().enumerate() {
                *counter += 1;
                let fresh = if (bits >> k) & 1 == 0 {
                    format!("o{}", counter)
                } else {
                    format!("v{}", counter)
                };
                mapping.insert(label.clone(), fresh);
            }
            result.push(term.renamed(&mapping));
        }
    }

    result
}
