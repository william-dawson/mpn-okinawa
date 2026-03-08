use okinawa_wick::step_1_expand::expand_general_labels;
use okinawa_wick::step_2_normal_order::normal_order_fermi_vacuum;
use okinawa_wick::step_3_filter::filter_fully_contracted;
use okinawa_wick::step_4_deltas::apply_deltas;
use okinawa_wick::step_4b_reclassify::reclassify_occ_repulsion;
use okinawa_wick::step_5_labels::canonicalize_labels_pdaggerq_style;
use okinawa_wick::step_6_cleanup::cancel_terms;
use okinawa_wick::step_6a_linked_cluster::filter_unlinked_diagrams;
use okinawa_wick::term::{is_general, is_occ, Delta, Operator, Tensor, Term};
use okinawa_wick::WickHelper;

#[test]
fn test_expand_general_labels_basic() {
    let t = Term::new(
        1.0,
        1,
        vec![Tensor::new("g", vec!["p", "q", "i", "a"])],
        vec![Operator::new("p", true), Operator::new("q", false)],
        vec![],
    );
    let out = expand_general_labels(&[t], &mut 0);

    assert_eq!(out.len(), 4);
    for term in &out {
        for tensor in &term.tensors {
            for idx in &tensor.indices {
                assert!(!is_general(idx));
            }
        }
        for op in &term.operators {
            assert!(!is_general(&op.label));
        }
    }
}

#[test]
fn test_normal_order_then_filter_fully_contracted() {
    let t = Term::new(
        1.0,
        1,
        vec![Tensor::new("g", vec!["i", "j", "i", "j"])],
        vec![Operator::new("i", true), Operator::new("j", false)],
        vec![],
    );
    let ordered = normal_order_fermi_vacuum(&[t]);
    let contracted = filter_fully_contracted(&ordered);

    assert_eq!(ordered.len(), 2);
    assert_eq!(contracted.len(), 1);
    assert!(contracted[0].operators.is_empty());
    assert_eq!(contracted[0].deltas.len(), 1);
}

#[test]
fn test_apply_deltas_eliminates_delta_and_unifies_label_pair() {
    let t = Term::new(
        1.0,
        1,
        vec![Tensor::new("g", vec!["o1", "v1", "o2", "v2"])],
        vec![Operator::new("o1", true)],
        vec![Delta::new("o1", "o2")],
    );
    let out = &apply_deltas(&[t])[0];

    assert!(out.deltas.is_empty());
    let labels_after: std::collections::HashSet<&String> = out
        .tensors
        .iter()
        .flat_map(|tensor| tensor.indices.iter())
        .filter(|idx| *idx == "o1" || *idx == "o2")
        .collect();
    assert_eq!(labels_after.len(), 1);
    assert!(out.operators[0].label == "o1" || out.operators[0].label == "o2");
}

#[test]
fn test_reclassify_occ_repulsion_into_eri() {
    let t = Term::new(
        1.0,
        1,
        vec![Tensor::from_strings(
            "occ_repulsion".to_string(),
            vec!["o1".to_string(), "o2".to_string()],
        )],
        vec![],
        vec![],
    );
    let out = &reclassify_occ_repulsion(&[t])[0];

    assert_eq!(out.tensors.len(), 1);
    let tensor = &out.tensors[0];
    assert_eq!(tensor.name, "g");
    assert_eq!(tensor.indices[0], "o1");
    assert_eq!(tensor.indices[2], "o2");
    assert_eq!(tensor.indices[1], tensor.indices[3]);
    assert!(is_occ(&tensor.indices[1]));
}

#[test]
fn test_canonicalize_labels_removes_fresh_prefixes() {
    let t = Term::new(
        1.0,
        1,
        vec![Tensor::new("g", vec!["o3", "v9", "o4", "v10"])],
        vec![],
        vec![],
    );
    let out = &canonicalize_labels_pdaggerq_style(&[t])[0];

    for idx in &out.tensors[0].indices {
        let is_fresh_occ =
            idx.starts_with('o') && idx.len() > 1 && idx[1..].chars().all(|c| c.is_ascii_digit());
        let is_fresh_virt =
            idx.starts_with('v') && idx.len() > 1 && idx[1..].chars().all(|c| c.is_ascii_digit());
        assert!(!is_fresh_occ);
        assert!(!is_fresh_virt);
    }
}

#[test]
fn test_cancel_terms_sums_and_cancels_identical_terms() {
    let t1 = Term::new(
        0.25,
        1,
        vec![Tensor::new("g", vec!["i", "j", "a", "b"])],
        vec![],
        vec![],
    );
    let t2 = Term::new(
        0.25,
        -1,
        vec![Tensor::new("g", vec!["i", "j", "a", "b"])],
        vec![],
        vec![],
    );
    let out = cancel_terms(&[t1, t2]);
    assert!(out.is_empty());
}

#[test]
fn test_linked_filter_behaves_reasonably() {
    let connected = Term::new(
        1.0,
        1,
        vec![
            Tensor::new("g", vec!["i", "j", "a", "b"]),
            Tensor::new("g", vec!["a", "b", "k", "l"]),
        ],
        vec![],
        vec![],
    );
    let disconnected = Term::new(
        1.0,
        1,
        vec![
            Tensor::new("g", vec!["i", "j", "i", "j"]),
            Tensor::new("g", vec!["k", "l", "k", "l"]),
        ],
        vec![],
        vec![],
    );
    let linked = filter_unlinked_diagrams(&[connected.clone(), disconnected]);
    assert_eq!(linked.len(), 1);
}

// ── Smoke tests matching Python test_mp_smoke.py ────────────────────────

#[test]
fn test_mp1_term_counts() {
    let mut wh_parity = WickHelper::new(false);
    wh_parity.add_operator_product(1.0, &["v"]);
    wh_parity.simplify();
    assert_eq!(wh_parity.strings().len(), 1);

    let mut wh_default = WickHelper::default();
    wh_default.add_operator_product(1.0, &["v"]);
    wh_default.simplify();
    assert_eq!(wh_default.strings().len(), 1);
}

#[test]
fn test_mp2_term_counts() {
    let mut wh_parity = WickHelper::new(false);
    wh_parity.add_operator_product(1.0, &["v", "v"]);
    wh_parity.simplify();
    assert_eq!(wh_parity.strings().len(), 2);

    let mut wh_default = WickHelper::default();
    wh_default.add_operator_product(1.0, &["v", "v"]);
    wh_default.simplify();
    assert_eq!(wh_default.strings().len(), 1);
}

#[test]
fn test_mp3_term_counts() {
    let mut wh_parity = WickHelper::new(false);
    wh_parity.add_operator_product(1.0, &["v", "v", "v"]);
    wh_parity.simplify();
    assert_eq!(wh_parity.strings().len(), 5);

    let mut wh_default = WickHelper::default();
    wh_default.add_operator_product(1.0, &["v", "v", "v"]);
    wh_default.simplify();
    assert_eq!(wh_default.strings().len(), 3);
}

#[test]
fn test_denominators_align_with_output_terms() {
    let mut wh = WickHelper::new(false);
    wh.add_operator_product(1.0, &["v", "v"]);
    wh.simplify();

    let rows = wh.strings();
    let denoms = wh.denominators();
    assert_eq!(rows.len(), denoms.len());
    assert!(!rows.is_empty());
    assert!(denoms.iter().all(|d| !d.is_empty()));
}
