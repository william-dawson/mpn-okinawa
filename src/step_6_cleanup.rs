//! Step 6: Cleanup / Simplification with ERI Canonicalization
//!
//! 1. Canonicalize each integral within each tensor (sort indices, apply sign).
//! 2. Collect terms by canonical form and sum their numerical factors.
//! 3. Drop any term whose total factor is (effectively) zero.

use std::collections::HashMap;

use crate::term::{is_occ, is_virt, Tensor, Term};

const EPSILON: f64 = 1e-12;

fn canonicalize_tensor(tensor: &Tensor) -> (Tensor, i32) {
    if tensor.indices.len() != 4 {
        return (tensor.clone(), 1);
    }

    let mut p = tensor.indices[0].clone();
    let mut q = tensor.indices[1].clone();
    let mut r = tensor.indices[2].clone();
    let mut s = tensor.indices[3].clone();
    let mut sign = 1;

    if p > q {
        std::mem::swap(&mut p, &mut q);
        sign *= -1;
    }
    if r > s {
        std::mem::swap(&mut r, &mut s);
        sign *= -1;
    }

    (
        Tensor::from_strings(tensor.name.clone(), vec![p, q, r, s]),
        sign,
    )
}

type TensorKey = (String, Vec<String>);
type TermCanonKey = Vec<TensorKey>;

fn normalize_labels(tensors: &[(String, Vec<String>)]) -> (TermCanonKey, i32) {
    let mut occ_names = "ijklmn".chars();
    let mut virt_names = "abcdef".chars();
    let mut mapping: HashMap<String, String> = HashMap::new();
    let mut sign = 1;

    let map_label = |lbl: &str,
                     mapping: &mut HashMap<String, String>,
                     occ_names: &mut std::str::Chars,
                     virt_names: &mut std::str::Chars|
     -> String {
        if let Some(mapped) = mapping.get(lbl) {
            return mapped.clone();
        }
        let new_label = if is_occ(lbl) {
            occ_names.next().unwrap().to_string()
        } else if is_virt(lbl) {
            virt_names.next().unwrap().to_string()
        } else {
            lbl.to_string()
        };
        mapping.insert(lbl.to_string(), new_label.clone());
        new_label
    };

    let mut normalized: Vec<(String, Vec<String>)> = Vec::new();
    for (name, idxs) in tensors {
        let mapped: Vec<String> = idxs
            .iter()
            .map(|x| map_label(x, &mut mapping, &mut occ_names, &mut virt_names))
            .collect();
        if mapped.len() == 4 {
            let t = Tensor::from_strings(name.clone(), mapped);
            let (ct, s) = canonicalize_tensor(&t);
            sign *= s;
            normalized.push((name.clone(), ct.indices));
        } else {
            normalized.push((name.clone(), mapped));
        }
    }
    normalized.sort();
    (normalized, sign)
}

fn term_key(term: &Term) -> (TermCanonKey, i32) {
    let mut occ_labels: Vec<String> = Vec::new();
    let mut virt_labels: Vec<String> = Vec::new();

    for t in &term.tensors {
        for idx in &t.indices {
            if is_occ(idx) && !occ_labels.contains(idx) {
                occ_labels.push(idx.clone());
            }
            if is_virt(idx) && !virt_labels.contains(idx) {
                virt_labels.push(idx.clone());
            }
        }
    }
    occ_labels.sort();
    virt_labels.sort();

    let mut best_key: Option<TermCanonKey> = None;
    let mut best_sign = 1;

    let occ_perms = permutations(&occ_labels);
    let virt_perms = permutations(&virt_labels);

    for occ_perm in &occ_perms {
        let mut occ_map: HashMap<String, String> = HashMap::new();
        for (orig, perm) in occ_labels.iter().zip(occ_perm.iter()) {
            occ_map.insert(orig.clone(), perm.clone());
        }

        for virt_perm in &virt_perms {
            let mut rename_map = occ_map.clone();
            for (orig, perm) in virt_labels.iter().zip(virt_perm.iter()) {
                rename_map.insert(orig.clone(), perm.clone());
            }

            let mut canon_tensors: Vec<(String, Vec<String>)> = Vec::new();
            let mut total_sign = 1;

            for tensor in &term.tensors {
                let renamed = tensor.renamed(&rename_map);
                let (canon_tensor, s) = canonicalize_tensor(&renamed);
                canon_tensors.push((canon_tensor.name.clone(), canon_tensor.indices.clone()));
                total_sign *= s;
            }

            canon_tensors.sort();
            let (key, s2) = normalize_labels(&canon_tensors);
            total_sign *= s2;

            let replace = match &best_key {
                None => true,
                Some(bk) => key < *bk || (key == *bk && total_sign > best_sign),
            };
            if replace {
                best_key = Some(key);
                best_sign = total_sign;
            }
        }
    }

    (best_key.unwrap_or_default(), best_sign)
}

fn permutations(items: &[String]) -> Vec<Vec<String>> {
    if items.is_empty() {
        return vec![vec![]];
    }
    let mut result = Vec::new();
    for i in 0..items.len() {
        let mut rest: Vec<String> = items.to_vec();
        let item = rest.remove(i);
        for mut perm in permutations(&rest) {
            perm.insert(0, item.clone());
            result.push(perm);
        }
    }
    result
}

pub fn cancel_terms(terms: &[Term]) -> Vec<Term> {
    let mut buckets: Vec<(TermCanonKey, f64)> = Vec::new();
    let mut key_index: HashMap<TermCanonKey, usize> = HashMap::new();

    for term in terms {
        let (key, canon_sign) = term_key(term);
        let effective = term.effective_factor() * canon_sign as f64;

        if let Some(&idx) = key_index.get(&key) {
            buckets[idx].1 += effective;
        } else {
            key_index.insert(key.clone(), buckets.len());
            buckets.push((key, effective));
        }
    }

    let mut result = Vec::new();
    for (key, total_factor) in &buckets {
        if total_factor.abs() > EPSILON {
            let tensors: Vec<Tensor> = key
                .iter()
                .map(|(name, indices)| Tensor::from_strings(name.clone(), indices.clone()))
                .collect();
            let sign = if *total_factor > 0.0 { 1 } else { -1 };
            result.push(Term::new(total_factor.abs(), sign, tensors, vec![], vec![]));
        }
    }

    result
}
