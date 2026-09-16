//! CLI: runs the Peptonizer2000 pipeline over a peptide-protein relationship TSV and writes a
//! `protein_id, probability` TSV, restoring the original protein names on the way out.

use std::error::Error;

const ALPHAS: &[f32] = &[0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99];
const BETAS: &[f32] = &[0.4, 0.5, 0.7, 0.9];
const PRIORS: &[f32] = &[0.3, 0.5, 0.7];
const PROTEINS_TO_RETURN: usize = 100;

fn main() {
    if peptonizer_analysis::has_help_flag() {
        println!("{}", usage());
        return;
    }

    if let Err(error) = run() {
        eprintln!("protein_inference: {error}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn Error>> {
    let arguments =
        peptonizer_analysis::parse_arguments("--peptide-proteins", "protein_inference_results.tsv")?;
    let (relationships, proteins_by_id) =
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
        PROTEINS_TO_RETURN,
        None,
    )?;
    let probabilities =
        peptonizer_analysis::restore_original_ids(result.probabilities, &proteins_by_id)?;
    peptonizer_analysis::write_results(arguments.output, "protein_id", probabilities)?;
    println!(
        "Selected parameter set: alpha={}, beta={}, prior={}",
        result.alpha, result.beta, result.prior
    );
    Ok(())
}

fn usage() -> &'static str {
    "Usage: protein_inference --peptide-proteins FILE --peptide-scores FILE --peptide-counts FILE [--output FILE]"
}