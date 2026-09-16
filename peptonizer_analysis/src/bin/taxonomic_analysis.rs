//! CLI: runs the Peptonizer2000 pipeline over a peptide-taxon relationship TSV, normalizing
//! taxon IDs to [`TAXONOMIC_RANK`] before weighing, and writes a `taxon_id, probability` TSV.

use std::error::Error;

const ALPHAS: &[f32] = &[0.8, 0.9, 0.99];
const BETAS: &[f32] = &[0.6, 0.7, 0.8, 0.9];
const PRIORS: &[f32] = &[0.3, 0.5];
const TAXA_TO_RETURN: usize = 100;
const TAXONOMIC_RANK: &str = "species";

fn main() {
    if peptonizer_analysis::has_help_flag() {
        println!("{}", usage());
        return;
    }

    if let Err(error) = run() {
        eprintln!("taxonomic_analysis: {error}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn Error>> {
    let arguments =
        peptonizer_analysis::parse_arguments("--peptide-taxa", "taxonomic_analysis_results.tsv")?;
    let relationships = peptonizer_analysis::read_relationships(arguments.relationships)?;
    let scores = peptonizer_analysis::read_scores(arguments.scores)?;
    let counts = peptonizer_analysis::read_counts(arguments.counts)?;
    let result = peptonizer_analysis::run_analysis(
        relationships,
        scores,
        counts,
        ALPHAS,
        BETAS,
        PRIORS,
        TAXA_TO_RETURN,
        Some(TAXONOMIC_RANK),
    )?;
    peptonizer_analysis::write_results(arguments.output, "taxon_id", result.probabilities)?;
    println!(
        "Selected parameter set: alpha={}, beta={}, prior={}",
        result.alpha, result.beta, result.prior
    );
    Ok(())
}

fn usage() -> &'static str {
    "Usage: taxonomic_analysis --peptide-taxa FILE --peptide-scores FILE --peptide-counts FILE [--output FILE]"
}