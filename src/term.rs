//! Core data structures for second-quantized operator algebra.
//!
//! Operator input convention (TRUE VACUUM):
//!   '+p'  -> a†_p  (creation operator)
//!   '-p'  -> a_p   (annihilation operator)
//!
//! Index label conventions:
//!   Occupied : i, j, k, l, m, n  OR  o1, o2, o3, ...
//!   Virtual  : a, b, c, d, e, f  OR  v1, v2, v3, ...
//!   General  : p, q, r, s, t, u, w  (must be expanded before use)

use std::collections::HashMap;
use std::fmt;

// ── index classification ────────────────────────────────────────────────────

const OCC_CHARS: &[char] = &['i', 'j', 'k', 'l', 'm', 'n'];
const VIRT_CHARS: &[char] = &['a', 'b', 'c', 'd', 'e', 'f'];

pub fn is_occ(label: &str) -> bool {
    let first = label.chars().next().unwrap();
    if OCC_CHARS.contains(&first) {
        return true;
    }
    first == 'o' && label.len() > 1 && label[1..].chars().all(|c| c.is_ascii_digit())
}

pub fn is_virt(label: &str) -> bool {
    let first = label.chars().next().unwrap();
    if VIRT_CHARS.contains(&first) {
        return true;
    }
    first == 'v' && label.len() > 1 && label[1..].chars().all(|c| c.is_ascii_digit())
}

pub fn is_general(label: &str) -> bool {
    !is_occ(label) && !is_virt(label)
}

// ── Operator ────────────────────────────────────────────────────────────────

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Operator {
    pub label: String,
    pub dagger: bool,
}

impl Operator {
    pub fn new(label: &str, dagger: bool) -> Self {
        Self {
            label: label.to_string(),
            dagger,
        }
    }

    /// Creator relative to Fermi vacuum?
    /// dagger_fermi = dagger XOR is_occ(label)
    pub fn dagger_fermi(&self) -> bool {
        if is_occ(&self.label) {
            return !self.dagger;
        }
        if is_virt(&self.label) {
            return self.dagger;
        }
        panic!(
            "dagger_fermi undefined for general label '{}'. Expand general labels first.",
            self.label
        );
    }

    pub fn renamed(&self, mapping: &HashMap<String, String>) -> Self {
        Self {
            label: mapping.get(&self.label).cloned().unwrap_or_else(|| self.label.clone()),
            dagger: self.dagger,
        }
    }
}

impl fmt::Display for Operator {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", if self.dagger { "+" } else { "-" }, self.label)
    }
}

impl fmt::Debug for Operator {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

// ── Tensor ──────────────────────────────────────────────────────────────────

#[derive(Clone, PartialEq, Eq)]
pub struct Tensor {
    pub name: String,
    pub indices: Vec<String>,
}

impl Tensor {
    pub fn new(name: &str, indices: Vec<&str>) -> Self {
        Self {
            name: name.to_string(),
            indices: indices.into_iter().map(|s| s.to_string()).collect(),
        }
    }

    pub fn from_strings(name: String, indices: Vec<String>) -> Self {
        Self { name, indices }
    }

    pub fn renamed(&self, mapping: &HashMap<String, String>) -> Self {
        Self {
            name: self.name.clone(),
            indices: self
                .indices
                .iter()
                .map(|i| mapping.get(i).cloned().unwrap_or_else(|| i.clone()))
                .collect(),
        }
    }

    pub fn to_string_pdaggerq(&self) -> String {
        if self.name == "g" {
            format!(
                "<{},{}||{},{}>",
                self.indices[0], self.indices[1], self.indices[2], self.indices[3]
            )
        } else if self.name == "f" {
            format!("f({},{})", self.indices[0], self.indices[1])
        } else {
            format!(
                "{}({})",
                self.name,
                self.indices.join(",")
            )
        }
    }
}

impl fmt::Display for Tensor {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_string_pdaggerq())
    }
}

impl fmt::Debug for Tensor {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

// ── Delta ───────────────────────────────────────────────────────────────────

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Delta {
    pub i: String,
    pub j: String,
}

impl Delta {
    pub fn new(i: &str, j: &str) -> Self {
        Self {
            i: i.to_string(),
            j: j.to_string(),
        }
    }

    pub fn renamed(&self, mapping: &HashMap<String, String>) -> Self {
        Self {
            i: mapping.get(&self.i).cloned().unwrap_or_else(|| self.i.clone()),
            j: mapping.get(&self.j).cloned().unwrap_or_else(|| self.j.clone()),
        }
    }
}

impl fmt::Display for Delta {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "d({},{})", self.i, self.j)
    }
}

impl fmt::Debug for Delta {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

// ── Term ────────────────────────────────────────────────────────────────────

#[derive(Clone)]
pub struct Term {
    pub factor: f64,
    pub sign: i32,
    pub tensors: Vec<Tensor>,
    pub operators: Vec<Operator>,
    pub deltas: Vec<Delta>,
}

impl Term {
    pub fn new(
        factor: f64,
        sign: i32,
        tensors: Vec<Tensor>,
        operators: Vec<Operator>,
        deltas: Vec<Delta>,
    ) -> Self {
        Self {
            factor,
            sign,
            tensors,
            operators,
            deltas,
        }
    }

    pub fn effective_factor(&self) -> f64 {
        self.sign as f64 * self.factor
    }

    pub fn renamed(&self, mapping: &HashMap<String, String>) -> Self {
        Self {
            factor: self.factor,
            sign: self.sign,
            tensors: self.tensors.iter().map(|t| t.renamed(mapping)).collect(),
            operators: self.operators.iter().map(|op| op.renamed(mapping)).collect(),
            deltas: self.deltas.iter().map(|d| d.renamed(mapping)).collect(),
        }
    }

    pub fn is_normal_ordered_fermi(&self) -> bool {
        let mut seen_annihilator = false;
        for op in &self.operators {
            let df = op.dagger_fermi();
            if !df {
                seen_annihilator = true;
            } else if seen_annihilator {
                return false;
            }
        }
        true
    }

    pub fn is_fully_contracted(&self) -> bool {
        self.operators.is_empty()
    }

    pub fn has_general_labels(&self) -> bool {
        for op in &self.operators {
            if is_general(&op.label) {
                return true;
            }
        }
        for t in &self.tensors {
            for idx in &t.indices {
                if is_general(idx) {
                    return true;
                }
            }
        }
        false
    }
}

impl fmt::Display for Term {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut parts = vec![format!("{:+.6}", self.effective_factor())];
        for t in &self.tensors {
            parts.push(t.to_string_pdaggerq());
        }
        for d in &self.deltas {
            parts.push(format!("{}", d));
        }
        for op in &self.operators {
            parts.push(format!("{}", op));
        }
        write!(f, "{}", parts.join(" * "))
    }
}

impl fmt::Debug for Term {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

// ── Parsing helpers ─────────────────────────────────────────────────────────

pub fn parse_tensor(s: &str) -> Tensor {
    let s = s.trim();
    if let Some(paren_start) = s.find('(') {
        let name = &s[..paren_start];
        let rest = &s[paren_start + 1..s.len() - 1];
        let indices: Vec<String> = rest.split(',').map(|x| x.trim().to_string()).collect();
        Tensor::from_strings(name.to_string(), indices)
    } else {
        panic!("Cannot parse tensor string: '{}'", s);
    }
}

pub fn parse_op(s: &str) -> Operator {
    let s = s.trim();
    let first = s.chars().next().unwrap();
    match first {
        '+' => Operator::new(&s[1..], true),
        '-' => Operator::new(&s[1..], false),
        _ => panic!("Operator string must start with '+' or '-', got: '{}'", s),
    }
}
