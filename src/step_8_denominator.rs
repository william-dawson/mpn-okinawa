//! Step 8: Extract energy denominators from final terms.

use std::collections::HashSet;

use crate::term::{is_occ, is_virt, Term};

pub fn extract_denominator(term: &Term) -> String {
    let mut seen = HashSet::new();
    let mut occ_indices = Vec::new();
    let mut virt_indices = Vec::new();

    for tensor in &term.tensors {
        for idx in &tensor.indices {
            if seen.insert(idx.clone()) {
                if is_occ(idx) {
                    occ_indices.push(idx.clone());
                } else if is_virt(idx) {
                    virt_indices.push(idx.clone());
                }
            }
        }
    }

    let mut parts = Vec::new();
    for idx in &occ_indices {
        parts.push(format!("e_{}", idx));
    }
    for idx in &virt_indices {
        parts.push(format!("-e_{}", idx));
    }

    if parts.is_empty() {
        return "1".to_string();
    }

    let mut result = parts[0].clone();
    for part in &parts[1..] {
        if part.starts_with('-') {
            result.push_str(&format!(" {}", part));
        } else {
            result.push_str(&format!(" + {}", part));
        }
    }

    result
}
