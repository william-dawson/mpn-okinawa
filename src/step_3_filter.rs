//! Step 3: Contraction Filter
//!
//! Keep only fully-contracted terms (no remaining operators).

use crate::term::Term;

pub fn filter_fully_contracted(terms: &[Term]) -> Vec<Term> {
    terms
        .iter()
        .filter(|t| t.is_fully_contracted())
        .cloned()
        .collect()
}
