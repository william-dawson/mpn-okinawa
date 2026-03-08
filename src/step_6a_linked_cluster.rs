//! Step 6a: Linked-Cluster Theorem Filter
//!
//! Remove unlinked (disconnected) diagram contributions.
//! A diagram is LINKED if all tensors are connected via shared indices.

use std::collections::HashSet;

use crate::term::Term;

fn tensors_share_index(t1: &crate::term::Tensor, t2: &crate::term::Tensor) -> bool {
    let s1: HashSet<&String> = t1.indices.iter().collect();
    t2.indices.iter().any(|idx| s1.contains(idx))
}

fn is_connected(term: &Term) -> bool {
    let tensors = &term.tensors;
    if tensors.len() <= 1 {
        return true;
    }

    let mut visited = HashSet::new();
    visited.insert(0usize);
    let mut stack = vec![0usize];

    while let Some(i) = stack.pop() {
        for j in 0..tensors.len() {
            if visited.contains(&j) {
                continue;
            }
            if tensors_share_index(&tensors[i], &tensors[j]) {
                visited.insert(j);
                stack.push(j);
            }
        }
    }

    visited.len() == tensors.len()
}

pub fn filter_unlinked_diagrams(terms: &[Term]) -> Vec<Term> {
    terms.iter().filter(|t| is_connected(t)).cloned().collect()
}
