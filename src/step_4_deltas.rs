//! Step 4: Delta Elimination
//!
//! Each Delta(i, j) means delta_{ij}. Substitute one label for the other
//! everywhere and remove the delta.

use std::collections::HashMap;

use crate::term::Term;

fn count_in_tensors_and_ops(term: &Term, label: &str) -> usize {
    let mut count = 0;
    for t in &term.tensors {
        count += t.indices.iter().filter(|idx| idx.as_str() == label).count();
    }
    for op in &term.operators {
        if op.label == label {
            count += 1;
        }
    }
    count
}

fn apply_one_delta(term: &Term) -> Term {
    if term.deltas.is_empty() {
        return term.clone();
    }

    let delta = &term.deltas[0];
    let remaining_deltas: Vec<_> = term.deltas[1..].to_vec();

    let i = &delta.i;
    let j = &delta.j;

    if i == j {
        return Term::new(
            term.factor,
            term.sign,
            term.tensors.clone(),
            term.operators.clone(),
            remaining_deltas,
        );
    }

    let (keep, replace) = if count_in_tensors_and_ops(term, i) > 0 {
        (j.clone(), i.clone())
    } else if count_in_tensors_and_ops(term, j) > 0 {
        (i.clone(), j.clone())
    } else {
        (j.clone(), i.clone())
    };

    let mut mapping = HashMap::new();
    mapping.insert(replace, keep);

    Term::new(
        term.factor,
        term.sign,
        term.tensors.iter().map(|t| t.renamed(&mapping)).collect(),
        term.operators.iter().map(|op| op.renamed(&mapping)).collect(),
        remaining_deltas.iter().map(|d| d.renamed(&mapping)).collect(),
    )
}

pub fn apply_deltas(terms: &[Term]) -> Vec<Term> {
    let mut result = Vec::new();
    for term in terms {
        let mut t = term.clone();
        while !t.deltas.is_empty() {
            t = apply_one_delta(&t);
        }
        result.push(t);
    }
    result
}
