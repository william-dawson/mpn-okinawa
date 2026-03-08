//! WickHelper: main entry point mirroring the pq_helper interface.

use crate::step_1_expand::expand_general_labels;
use crate::step_2_normal_order::normal_order_fermi_vacuum;
use crate::step_3_filter::filter_fully_contracted;
use crate::step_4_deltas::apply_deltas;
use crate::step_4b_reclassify::reclassify_occ_repulsion;
use crate::step_5_labels::canonicalize_labels_pdaggerq_style;
use crate::step_6_cleanup::cancel_terms;
use crate::step_6a_linked_cluster::filter_unlinked_diagrams;
use crate::step_7_output::format_strings;
use crate::step_8_denominator::extract_denominator;
use crate::term::{parse_op, parse_tensor, Operator, Tensor, Term};

pub struct WickHelper {
    input_terms: Vec<Term>,
    result: Vec<Term>,
    filter_unlinked: bool,
}

impl WickHelper {
    pub fn new(filter_unlinked: bool) -> Self {
        Self {
            input_terms: Vec::new(),
            result: Vec::new(),
            filter_unlinked,
        }
    }

    pub fn add_term(&mut self, factor: f64, tensors: &[&str], ops: &[&str]) {
        let parsed_tensors: Vec<Tensor> = tensors.iter().map(|t| parse_tensor(t)).collect();
        let parsed_ops: Vec<Operator> = ops.iter().map(|o| parse_op(o)).collect();
        self.input_terms
            .push(Term::new(factor, 1, parsed_tensors, parsed_ops, vec![]));
    }

    pub fn add_operator_product(&mut self, factor: f64, ops: &[&str]) {
        for op in ops {
            if op.to_lowercase() != "v" {
                panic!("add_operator_product currently supports only 'v' inputs");
            }
        }

        let mut branches: Vec<(f64, Vec<Tensor>, Vec<Operator>)> =
            vec![(factor, vec![], vec![])];
        let mut next_idx = 0usize;

        for _ in ops {
            let mut new_branches = Vec::new();
            for (coef, tensors, op_list) in &branches {
                let p1 = format!("p{}", next_idx);
                let p2 = format!("p{}", next_idx + 1);
                let p3 = format!("p{}", next_idx + 2);
                let p4 = format!("p{}", next_idx + 3);

                // j1 piece
                {
                    let mut new_tensors = tensors.clone();
                    new_tensors.push(Tensor::from_strings(
                        "occ_repulsion".to_string(),
                        vec![p1.clone(), p2.clone()],
                    ));
                    let mut new_ops = op_list.clone();
                    new_ops.push(Operator::new(&p1, true));
                    new_ops.push(Operator::new(&p2, false));
                    new_branches.push((-coef, new_tensors, new_ops));
                }

                // j2 piece
                {
                    let mut new_tensors = tensors.clone();
                    new_tensors.push(Tensor::from_strings(
                        "g".to_string(),
                        vec![p1.clone(), p2.clone(), p4.clone(), p3.clone()],
                    ));
                    let mut new_ops = op_list.clone();
                    new_ops.push(Operator::new(&p1, true));
                    new_ops.push(Operator::new(&p2, true));
                    new_ops.push(Operator::new(&p3, false));
                    new_ops.push(Operator::new(&p4, false));
                    new_branches.push((0.25 * coef, new_tensors, new_ops));
                }
            }
            branches = new_branches;
            next_idx += 4;
        }

        for (coef, tensors, op_list) in branches {
            let sign = if coef >= 0.0 { 1 } else { -1 };
            self.input_terms
                .push(Term::new(coef.abs(), sign, tensors, op_list, vec![]));
        }
    }

    pub fn simplify(&mut self) {
        let mut terms = self.input_terms.clone();
        let mut counter = 0usize;

        // 1. Expand general labels
        terms = expand_general_labels(&terms, &mut counter);

        // 2. Normal ordering (Wick's theorem)
        terms = normal_order_fermi_vacuum(&terms);

        // 3. Keep only fully contracted terms
        terms = filter_fully_contracted(&terms);

        // 4. Apply delta functions
        terms = apply_deltas(&terms);

        // 4b. Reclassify occ-repulsion intermediates
        terms = reclassify_occ_repulsion(&terms);

        // 5. Canonicalize index labels
        terms = canonicalize_labels_pdaggerq_style(&terms);

        // 6. Consolidate duplicate terms
        terms = cancel_terms(&terms);

        // 6a. Optional: filter unlinked diagrams
        if self.filter_unlinked {
            terms = filter_unlinked_diagrams(&terms);
        }

        self.result = terms;
    }

    pub fn strings(&self) -> Vec<Vec<String>> {
        format_strings(&self.result)
    }

    pub fn denominators(&self) -> Vec<String> {
        self.result.iter().map(|t| extract_denominator(t)).collect()
    }

    pub fn terms_with_denominators(&self) -> Vec<(Vec<String>, String)> {
        let strings = self.strings();
        let denoms = self.denominators();
        strings.into_iter().zip(denoms).collect()
    }

    pub fn result_terms(&self) -> &[Term] {
        &self.result
    }

    pub fn take_input_terms(&mut self) -> Vec<Term> {
        std::mem::take(&mut self.input_terms)
    }

    pub fn clear(&mut self) {
        self.input_terms.clear();
        self.result.clear();
    }
}

impl Default for WickHelper {
    fn default() -> Self {
        Self::new(true)
    }
}
