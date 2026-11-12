//! CLI: runs the Peptonizer2000 pipeline over a peptide-function relationship TSV and writes a
//! `function_id, probability` TSV. Function IDs may be arbitrary strings (e.g. GO or EC terms),
//! restored on the way out; the input file must have no header row.

use std::error::Error;

const ALPHAS: &[f32] = &[0.8, 0.9, 0.99];
const BETAS: &[f32] = &[0.6, 0.7, 0.8, 0.9];
const PRIORS: &[f32] = &[0.3, 0.5];
const DEFAULT_FUNCTIONS_TO_RETURN: usize = 100;

fn main() {
    if peptonizer_analysis::has_help_flag() {
        println!("{}", usage());
        return;
    }

    if let Err(error) = run() {
        eprintln!("functional_analysis: {error}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn Error>> {
    let arguments = peptonizer_analysis::parse_arguments(
        "--peptide-functions",
        "--functions-to-return",
        "functional_analysis_results.tsv",
        DEFAULT_FUNCTIONS_TO_RETURN,
    )?;
    let (relationships, functions_by_id) =
        peptonizer_analysis::read_relationships_with_string_ids(arguments.relationships)?;
    let scores = peptonizer_analysis::read_scores(arguments.scores)?;
    let counts = peptonizer_analysis::read_counts(arguments.counts)?;
    let result = peptonizer_analysis::run_analysis(
        relationships,
        scores,
        counts,
        ALPHAS,
        BETAS,
        PRIORS,
        arguments.results_to_return,
        arguments.max_peptides,
        None,
    )?;
    let probabilities =
        peptonizer_analysis::restore_original_ids(result.probabilities, &functions_by_id)?;
    peptonizer_analysis::write_results(arguments.output, "function_id", probabilities)?;
    println!(
        "Selected parameter set: alpha={}, beta={}, prior={} (functions-to-return={}, max-peptides={})",
        result.alpha, result.beta, result.prior, arguments.results_to_return, arguments.max_peptides
    );
    Ok(())
}

fn usage() -> String {
    format!(
        "Usage: functional_analysis --peptide-functions FILE --peptide-scores FILE --peptide-counts FILE [--output FILE] [--functions-to-return N (default {})] [--max-peptides N (default {})]",
        DEFAULT_FUNCTIONS_TO_RETURN,
        peptonizer_analysis::DEFAULT_MAX_PEPTIDES,
    )
}