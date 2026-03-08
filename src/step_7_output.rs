//! Step 7: Output Formatting
//!
//! Convert fully-contracted Terms into pdaggerq strings() format.

use crate::term::Term;

fn fmt_factor(value: f64) -> String {
    let sign = if value >= 0.0 { '+' } else { '-' };
    let abs_val = value.abs();
    let scaled = 10.0 * abs_val;

    let s = format!("{:.6}", scaled);

    let mut precision: usize = 0;
    let mut decimal_point_found = false;
    let mut last_digit = ' ';
    let mut repeat_count: usize = 0;
    let mut is_repeated;

    for digit in s.chars() {
        is_repeated = digit == last_digit;
        last_digit = digit;

        if digit == '.' {
            decimal_point_found = true;
        } else if decimal_point_found && is_repeated {
            repeat_count += 1;
            if repeat_count >= 12 {
                break;
            }
        }

        if !is_repeated {
            repeat_count = 0;
        }

        if decimal_point_found {
            precision += 1;
        }
    }

    if precision >= repeat_count && last_digit == '0' {
        precision -= repeat_count;
    }

    if precision < 2 {
        precision = 2;
    }

    format!("{}{:.*}", sign, precision, abs_val)
}

pub fn format_strings(terms: &[Term]) -> Vec<Vec<String>> {
    let mut result = Vec::new();
    for term in terms {
        if !term.is_fully_contracted() {
            continue;
        }
        let mut row = vec![fmt_factor(term.effective_factor())];
        for t in &term.tensors {
            row.push(t.to_string_pdaggerq());
        }
        result.push(row);
    }
    result
}
