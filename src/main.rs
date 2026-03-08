//! Wick's Theorem Engine CLI.

use clap::Parser;

use okinawa_wick::step_1_expand::expand_general_labels;
use okinawa_wick::step_2_normal_order::normal_order_fermi_vacuum;
use okinawa_wick::step_3_filter::filter_fully_contracted;
use okinawa_wick::step_4_deltas::apply_deltas;
use okinawa_wick::step_4b_reclassify::reclassify_occ_repulsion;
use okinawa_wick::step_5_labels::canonicalize_labels_pdaggerq_style;
use okinawa_wick::step_6_cleanup::cancel_terms;
use okinawa_wick::step_6a_linked_cluster::filter_unlinked_diagrams;
use okinawa_wick::step_7_output::format_strings;
use okinawa_wick::step_8_denominator::extract_denominator;
use okinawa_wick::term::Term;
use okinawa_wick::WickHelper;

#[derive(Parser)]
#[command(
    name = "wick",
    about = "Wick's Theorem Engine CLI",
    after_help = "Examples:\n  wick --mp1\n  wick --mp2 --verbose\n  wick --mp 5 --verbose\n  wick --mp12\n  wick --mp123"
)]
struct Cli {
    #[arg(long, help = "Run MP1 example")]
    mp1: bool,

    #[arg(long, help = "Run MP2 example")]
    mp2: bool,

    #[arg(long, help = "Run MP3 example")]
    mp3: bool,

    #[arg(long, help = "Run MP4 example (slow, 65K+ terms)")]
    mp4: bool,

    #[arg(long, help = "Run generic MPn (e.g. --mp 5)")]
    mp: Option<usize>,

    #[arg(long, help = "Run MP1+MP2 combined example")]
    mp12: bool,

    #[arg(long, help = "Run MP1+MP2+MP3 combined example")]
    mp123: bool,

    #[arg(short, long, help = "Show detailed output for each pipeline step")]
    verbose: bool,
}

fn print_section(title: &str) {
    println!("\n{}", "=".repeat(72));
    println!(" {}", title);
    println!("{}\n", "=".repeat(72));
}

fn print_terms(terms: &[Term], limit: Option<usize>) {
    if terms.is_empty() {
        println!("  (empty)");
        return;
    }

    let shown = match limit {
        Some(n) => &terms[..n.min(terms.len())],
        None => terms,
    };

    for (i, term) in shown.iter().enumerate() {
        let factor_str = if term.sign > 0 {
            format!("{:+.6}", term.factor)
        } else {
            format!("{:+.6}", -term.factor)
        };
        let ops_str: Vec<String> = term
            .operators
            .iter()
            .map(|op| {
                format!(
                    "{}{}",
                    op.label,
                    if op.dagger { "\u{2020}" } else { "" }
                )
            })
            .collect();
        let deltas_str: Vec<String> = term
            .deltas
            .iter()
            .map(|d| format!("\u{03B4}({},{})", d.i, d.j))
            .collect();
        let tensors_str: Vec<String> = term.tensors.iter().map(|t| t.to_string_pdaggerq()).collect();

        println!("  [{}] factor={}, sign={:+}", i, factor_str, term.sign);
        if !tensors_str.is_empty() {
            println!("       tensors: [{}]", tensors_str.join(", "));
        }
        if !ops_str.is_empty() {
            println!("       ops: [{}]", ops_str.join(", "));
        }
        if !deltas_str.is_empty() {
            println!("       deltas: [{}]", deltas_str.join(", "));
        }
    }

    if let Some(n) = limit {
        if terms.len() > n {
            println!("  ... and {} more terms", terms.len() - n);
        }
    }
    println!();
}

fn print_strings_with_denominators(
    strings_output: &[Vec<String>],
    denominators: &[String],
    limit: Option<usize>,
) {
    if strings_output.is_empty() {
        println!("  (empty)");
        return;
    }

    let n = match limit {
        Some(n) => n.min(strings_output.len()),
        None => strings_output.len(),
    };

    for (i, (row, denom)) in strings_output[..n]
        .iter()
        .zip(&denominators[..n])
        .enumerate()
    {
        let factor = &row[0];
        let tensors = if row.len() > 1 {
            row[1..].join(" * ")
        } else {
            "(scalar)".to_string()
        };
        println!("  [{:>2}] {:>8} * {}", i, factor, tensors);
        println!("         / ({})", denom);
    }

    if let Some(lim) = limit {
        if strings_output.len() > lim {
            println!("  ... and {} more terms", strings_output.len() - lim);
        }
    }
    println!();
}

fn build_input_terms(order: usize) -> Vec<Term> {
    let mut wh = WickHelper::new(true);
    let v_ops: Vec<&str> = (0..order).map(|_| "v").collect();
    wh.add_operator_product(1.0, &v_ops);
    // We need access to input terms; simplify builds them internally.
    // Instead, reconstruct the terms the same way add_operator_product does.
    // Actually, let's expose input_terms via a method.
    // For now, we'll just call simplify and show final results in verbose mode.
    // Let's use a different approach - rebuild the pipeline manually.
    drop(wh);

    let mut helper = WickHelper::new(true);
    helper.add_operator_product(1.0, &v_ops);
    helper.take_input_terms()
}

fn run_verbose_pipeline(terms: Vec<Term>) -> (Vec<Vec<String>>, Vec<String>) {
    let mut terms = terms;

    print_section("Step 1: General Label Expansion");
    let mut counter = 0;
    terms = expand_general_labels(&terms, &mut counter);
    println!("Output: {} terms", terms.len());
    print_terms(&terms, Some(3));

    print_section("Step 2: Normal Ordering (Wick's Theorem)");
    terms = normal_order_fermi_vacuum(&terms);
    println!("Output: {} terms", terms.len());
    print_terms(&terms, Some(3));

    print_section("Step 3: Filter Fully Contracted");
    terms = filter_fully_contracted(&terms);
    println!("Output: {} terms", terms.len());
    print_terms(&terms, Some(3));

    print_section("Step 4: Delta Elimination");
    terms = apply_deltas(&terms);
    println!("Output: {} terms", terms.len());
    print_terms(&terms, Some(3));

    print_section("Step 4b: Reclassify occ_repulsion");
    terms = reclassify_occ_repulsion(&terms);
    println!("Output: {} terms", terms.len());
    print_terms(&terms, Some(3));

    print_section("Step 5: Label Canonicalization");
    terms = canonicalize_labels_pdaggerq_style(&terms);
    println!("Output: {} terms", terms.len());
    print_terms(&terms, Some(3));

    print_section("Step 6: Consolidation");
    terms = cancel_terms(&terms);
    println!("Output: {} terms", terms.len());
    print_terms(&terms, Some(3));

    print_section("Step 6a: Linked-Diagram Filter");
    terms = filter_unlinked_diagrams(&terms);
    println!("Output: {} terms", terms.len());
    print_terms(&terms, Some(3));

    print_section("Step 7: Format Strings");
    let result = format_strings(&terms);
    println!("Output: {} terms", result.len());

    print_section("Step 8: Energy Denominator Extraction");
    let denoms: Vec<String> = terms.iter().map(|t| extract_denominator(t)).collect();
    println!("Output: {} denominators", denoms.len());
    print_strings_with_denominators(&result, &denoms, None);

    (result, denoms)
}

fn run_mpn(order: usize, verbose: bool) {
    if order < 1 {
        eprintln!("MP order must be >= 1");
        std::process::exit(1);
    }

    print_section(&format!("Running MP{}", order));
    println!("Input:");
    println!("  operator product = ['v'] * {}", order);

    print_section("Running Pipeline");
    if verbose {
        let terms = build_input_terms(order);
        run_verbose_pipeline(terms);
    } else {
        let mut wh = WickHelper::new(true);
        let v_ops: Vec<&str> = (0..order).map(|_| "v").collect();
        wh.add_operator_product(1.0, &v_ops);
        wh.simplify();
        let result = wh.strings();
        let denoms = wh.denominators();
        println!("Our output ({} terms) with energy denominators:", result.len());
        print_strings_with_denominators(&result, &denoms, None);
    }
}

fn run_combined(orders: &[usize], label: &str, verbose: bool) {
    print_section(&format!("Running {}", label));

    let mut wh = WickHelper::new(true);
    for &n in orders {
        let v_ops: Vec<&str> = (0..n).map(|_| "v").collect();
        wh.add_operator_product(1.0, &v_ops);
    }

    print_section("Running Pipeline");
    if verbose {
        let terms = wh.take_input_terms();
        let (result, denoms) = run_verbose_pipeline(terms);
        let _ = (result, denoms);
    } else {
        wh.simplify();
        let mut result = wh.strings();
        result.sort();
        let denoms = wh.denominators();
        println!(
            "Our combined output ({} terms) with energy denominators:",
            result.len()
        );
        print_strings_with_denominators(&result, &denoms, None);
    }
}

fn main() {
    let cli = Cli::parse();

    let has_command = cli.mp1
        || cli.mp2
        || cli.mp3
        || cli.mp4
        || cli.mp12
        || cli.mp123
        || cli.mp.is_some();

    if !has_command {
        use clap::CommandFactory;
        Cli::command().print_help().unwrap();
        println!();
        return;
    }

    if let Some(n) = cli.mp {
        run_mpn(n, cli.verbose);
    } else if cli.mp1 {
        run_mpn(1, cli.verbose);
    } else if cli.mp2 {
        run_mpn(2, cli.verbose);
    } else if cli.mp12 {
        run_combined(&[1, 2], "MP1+MP2 (Combined)", cli.verbose);
    } else if cli.mp123 {
        run_combined(&[1, 2, 3], "MP1+MP2+MP3 (Combined)", cli.verbose);
    } else if cli.mp3 {
        run_mpn(3, cli.verbose);
    } else if cli.mp4 {
        run_mpn(4, cli.verbose);
    }
}
