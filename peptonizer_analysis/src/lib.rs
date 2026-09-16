//! Native, standalone entry point for the Peptonizer2000 belief-propagation pipeline.
//!
//! Unlike the Snakemake workflow and the browser/WASM frontend, this crate runs the pipeline
//! directly from tab-separated peptide relationship/score/count files with no Python involved.
//! `protein_inference` and `functional_analysis` make no Unipept queries: protein and function
//! IDs are used as-is. `taxonomic_analysis` does query Unipept — it normalizes taxon IDs to
//! [`taxonomic_analysis`'s configured rank](../bin/taxonomic_analysis.rs) before weighing, the
//! same way the Snakemake workflow and browser frontend do, so it requires network access to
//! `api.unipept.ugent.be`. This crate reuses `peptonizer_rust`'s algorithm code (`weight_effects`,
//! `factor_graph`, `effects_clustering`, `zero_lookahead_belief_propagation`,
//! `analyse_grid_search`) as a normal path dependency — see that crate's `Cargo.toml` for the
//! `rlib` build target and `lib.rs` for the `pub mod` declarations that make this possible.
//!
//! The `taxonomic_analysis`, `protein_inference`, and `functional_analysis` binaries in `src/bin/`
//! are thin wrappers around [`run_analysis`] that only differ in their CLI flag name, default
//! output filename, and grid-search parameter ranges.

use std::collections::HashMap;
use std::error::Error;
use std::fs;
use std::future::Future;
use std::path::Path;
use std::sync::OnceLock;

use csv::{ReaderBuilder, WriterBuilder};
use peptonizer_rust::{analyse_grid_search, effects_clustering, factor_graph, weight_effects, zero_lookahead_belief_propagation};
use tokio::runtime::{Builder, Runtime};

/// Parsed command-line arguments shared by all three analysis binaries.
pub struct AnalysisArguments {
    /// Path to the peptide-relationship TSV (peptide → taxon/protein/function ID).
    pub relationships: String,
    /// Path to the peptide-score TSV (peptide → search-engine score).
    pub scores: String,
    /// Path to the peptide-count TSV (peptide → PSM count).
    pub counts: String,
    /// Path the results TSV should be written to.
    pub output: String,
}

/// Outcome of running the grid search: the winning posterior probabilities together with the
/// `(alpha, beta, prior)` parameter set that produced them.
pub struct AnalysisResult {
    /// Posterior probability per relationship ID, keyed by that ID's string representation.
    pub probabilities: HashMap<String, f32>,
    pub alpha: f32,
    pub beta: f32,
    pub prior: f32,
}

/// Returns `true` if `--help` or `-h` was passed on the command line.
pub fn has_help_flag() -> bool {
    has_help_flag_in(std::env::args())
}

fn has_help_flag_in(mut args: impl Iterator<Item = String>) -> bool {
    args.any(|argument| argument == "--help" || argument == "-h")
}

/// Parses the process's command-line arguments into an [`AnalysisArguments`].
///
/// `relationship_flag` is the binary-specific flag that carries the relationship TSV path (e.g.
/// `--peptide-taxa`, `--peptide-proteins`, `--peptide-functions`); `default_output` is used when
/// `--output` is not given. Every recognized flag must be followed by a value; unrecognized flags
/// are rejected.
///
/// # Errors
/// Returns an error if a flag is missing its value, an unrecognized flag is passed, or
/// `relationship_flag`, `--peptide-scores`, or `--peptide-counts` is not provided.
pub fn parse_arguments(
    relationship_flag: &str,
    default_output: &str,
) -> Result<AnalysisArguments, Box<dyn Error>> {
    parse_arguments_from(std::env::args(), relationship_flag, default_output)
}

fn parse_arguments_from(
    args: impl Iterator<Item = String>,
    relationship_flag: &str,
    default_output: &str,
) -> Result<AnalysisArguments, Box<dyn Error>> {
    let mut relationships = None;
    let mut scores = None;
    let mut counts = None;
    let mut output = default_output.to_owned();
    let mut arguments = args.skip(1);

    while let Some(argument) = arguments.next() {
        let value = arguments
            .next()
            .ok_or_else(|| format!("Missing value for {argument}"))?;
        match argument.as_str() {
            "--peptide-scores" => scores = Some(value),
            "--peptide-counts" => counts = Some(value),
            "--output" => output = value,
            _ if argument == relationship_flag => relationships = Some(value),
            _ => return Err(format!("Unknown argument: {argument}").into()),
        }
    }

    Ok(AnalysisArguments {
        relationships: relationships.ok_or_else(|| format!("Missing {relationship_flag}"))?,
        scores: scores.ok_or("Missing --peptide-scores")?,
        counts: counts.ok_or("Missing --peptide-counts")?,
        output,
    })
}

/// Reads a `peptide<TAB>relationship_id` TSV into a peptide → relationship-IDs map.
///
/// The ID column must be a numeric ID (e.g. an NCBI taxon ID) — use
/// [`read_relationships_with_string_ids`] instead for arbitrary string IDs (protein names,
/// GO/EC-style functional annotation IDs, ...).
///
/// A file may optionally start with a header row: if the second column of the first row does not
/// parse as a `usize`, that row is silently skipped. Every other row must parse, and a peptide can
/// appear on multiple rows to be associated with multiple IDs.
///
/// # Errors
/// Returns an error if a row is missing a column, or if a non-header row's ID column does not
/// parse as a `usize`.
pub fn read_relationships(
    path: impl AsRef<Path>,
) -> Result<HashMap<String, Vec<usize>>, Box<dyn Error>> {
    let mut relationships = HashMap::new();
    let mut reader = tsv_reader(path)?;
    reader.set_headers(csv::StringRecord::new());
    for (index, record) in reader.records().enumerate() {
        let record = record?;
        let peptide = record.get(0).ok_or("Missing peptide column")?.to_owned();
        let relationship = match record.get(1).ok_or("Missing ID column")?.parse::<usize>() {
            Ok(relationship) => relationship,
            Err(_) if index == 0 => continue,
            Err(error) => return Err(error.into()),
        };
        relationships
            .entry(peptide)
            .or_insert_with(Vec::new)
            .push(relationship);
    }
    Ok(relationships)
}

/// Reads a `peptide<TAB>id` TSV into a peptide → internal-numeric-ID map, along with the reverse
/// mapping from internal ID back to the original ID string. Used for relationship IDs that are
/// arbitrary strings rather than numbers — protein names (`protein_inference`) and functional
/// annotation IDs such as GO/EC terms (`functional_analysis`).
///
/// Internal IDs are assigned sequentially in first-seen order so they can be used as compact node
/// IDs in the belief-propagation graph; [`restore_original_ids`] reverses this mapping on the
/// final results. Unlike [`read_relationships`]/[`read_scores`]/[`read_counts`], this function does
/// **not** tolerate an optional header row — these IDs are arbitrary strings, so a header cannot
/// be distinguished from real data by a failed parse, and would otherwise be silently ingested as
/// a bogus entry. The input file must not have a header row.
///
/// # Errors
/// Returns an error if a row is missing its peptide or ID column.
pub fn read_relationships_with_string_ids(
    path: impl AsRef<Path>,
) -> Result<(HashMap<String, Vec<usize>>, HashMap<usize, String>), Box<dyn Error>> {
    let mut relationships = HashMap::new();
    let mut interned_ids = HashMap::new();
    let mut ids_by_index = HashMap::new();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(path)?;

    for record in reader.records() {
        let record = record?;
        let peptide = record.get(0).ok_or("Missing peptide column")?.to_owned();
        let id = record.get(1).ok_or("Missing ID column")?.to_owned();
        let next_index = interned_ids.len();
        let index = *interned_ids.entry(id.clone()).or_insert(next_index);

        ids_by_index.entry(index).or_insert(id);
        relationships
            .entry(peptide)
            .or_insert_with(Vec::new)
            .push(index);
    }

    Ok((relationships, ids_by_index))
}

/// Replaces the internal numeric IDs used as [`AnalysisResult::probabilities`] keys with the
/// original ID strings, using the reverse mapping produced by
/// [`read_relationships_with_string_ids`].
///
/// # Errors
/// Returns an error if a key does not parse as a `usize`, or if it has no corresponding entry in
/// `ids_by_index`.
pub fn restore_original_ids(
    results: HashMap<String, f32>,
    ids_by_index: &HashMap<usize, String>,
) -> Result<HashMap<String, f32>, Box<dyn Error>> {
    results
        .into_iter()
        .map(|(index, probability)| {
            let index = index.parse::<usize>()?;
            let id = ids_by_index
                .get(&index)
                .ok_or_else(|| format!("No ID found for internal index {index}"))?;
            Ok((id.clone(), probability))
        })
        .collect()
}

/// Reads a `peptide<TAB>score` TSV into a peptide → score map.
///
/// Tolerates an optional header row using the same rule as [`read_relationships`]: if the first
/// row's score column does not parse as an `f32`, that row is skipped.
///
/// # Errors
/// Returns an error if a row is missing a column, or if a non-header row's score column does not
/// parse as an `f32`.
pub fn read_scores(path: impl AsRef<Path>) -> Result<HashMap<String, f32>, Box<dyn Error>> {
    let mut scores = HashMap::new();
    let mut reader = tsv_reader(path)?;
    reader.set_headers(csv::StringRecord::new());
    for (index, record) in reader.records().enumerate() {
        let record = record?;
        let peptide = record.get(0).ok_or("Missing peptide column")?.to_owned();
        let score = match record.get(1).ok_or("Missing score column")?.parse::<f32>() {
            Ok(score) => score,
            Err(_) if index == 0 => continue,
            Err(error) => return Err(error.into()),
        };
        scores.insert(peptide, score);
    }
    Ok(scores)
}

/// Reads a `peptide<TAB>count` TSV into a peptide → PSM-count map.
///
/// Tolerates an optional header row using the same rule as [`read_relationships`]: if the first
/// row's count column does not parse as a `usize`, that row is skipped.
///
/// # Errors
/// Returns an error if a row is missing a column, or if a non-header row's count column does not
/// parse as a `usize`.
pub fn read_counts(path: impl AsRef<Path>) -> Result<HashMap<String, usize>, Box<dyn Error>> {
    let mut counts = HashMap::new();
    let mut reader = tsv_reader(path)?;
    reader.set_headers(csv::StringRecord::new());
    for (index, record) in reader.records().enumerate() {
        let record = record?;
        let peptide = record.get(0).ok_or("Missing peptide column")?.to_owned();
        let count = match record.get(1).ok_or("Missing count column")?.parse::<usize>() {
            Ok(count) => count,
            Err(_) if index == 0 => continue,
            Err(error) => return Err(error.into()),
        };
        counts.insert(peptide, count);
    }
    Ok(counts)
}

/// Runs the full pipeline — effect weighing, factor graph construction, a belief-propagation grid
/// search over every `(alpha, beta, prior)` combination in `alphas`/`betas`/`priors`, and
/// goodness-based selection of the best-scoring parameter set.
///
/// `effects_rank` requests Unipept-rank normalization of relationship IDs before weighing (used
/// for taxonomic analysis); pass `None` to use relationship IDs as-is (protein inference, and
/// functional analysis, have no such rank).
///
/// # Errors
/// Returns an error if the inputs fail [`validate_inputs`], if any pipeline stage fails, or if
/// `alphas`, `betas`, or `priors` is empty (no parameter set to select as "best").
#[allow(clippy::too_many_arguments)]
pub fn run_analysis(
    relationships: HashMap<String, Vec<usize>>,
    scores: HashMap<String, f32>,
    counts: HashMap<String, usize>,
    alphas: &[f32],
    betas: &[f32],
    priors: &[f32],
    results_to_return: usize,
    effects_rank: Option<&str>,
) -> Result<AnalysisResult, Box<dyn Error>> {
    validate_inputs(&relationships, &scores, &counts)?;

    println!("Preparing peptide relationships and weights...");
    let (sequence_scores_csv, relationship_weights_csv) = block_on(weight_effects::perform_effects_weighing_typed(
        relationships,
        scores,
        counts,
        results_to_return,
        effects_rank.map(str::to_owned),
    ))?;

    println!("Building the factor graph...");
    let factor_graph = factor_graph::generate_graph(sequence_scores_csv.clone())?;

    println!("Clustering related IDs for parameter selection...");
    let cluster_heads_csv =
        effects_clustering::cluster_effects(sequence_scores_csv, relationship_weights_csv, 0.9)?;

    let mut best_goodness = f64::NEG_INFINITY;
    let mut best_result = None;
    let mut best_parameters = None;
    let parameter_set_count = alphas.len() * betas.len() * priors.len();
    let mut parameter_set_index = 0;
    for &alpha in alphas {
        for &beta in betas {
            for &prior in priors {
                parameter_set_index += 1;
                println!(
                    "Running belief propagation for parameter set {parameter_set_index}/{parameter_set_count} (alpha={alpha}, beta={beta}, prior={prior})..."
                );
                let result_json = zero_lookahead_belief_propagation::run_belief_propagation(
                    &factor_graph,
                    alpha,
                    beta,
                    true,
                    prior,
                    None,
                    None,
                )?;
                println!("Scoring parameter set {parameter_set_index}/{parameter_set_count}...");
                let goodness = analyse_grid_search::compute_goodness(&cluster_heads_csv, &result_json)?;

                if goodness > best_goodness {
                    best_goodness = goodness;
                    best_result = Some(serde_json::from_str(&result_json)?);
                    best_parameters = Some((alpha, beta, prior));
                }
            }
        }
    }

    println!("Selecting the best-scoring parameter set...");
    let probabilities = best_result.ok_or("No parameter sets were provided")?;
    let (alpha, beta, prior) = best_parameters.ok_or("No parameter sets were provided")?;
    Ok(AnalysisResult {
        probabilities,
        alpha,
        beta,
        prior,
    })
}

/// Writes `results` to `path` as a TSV with a `[id_header, "probability"]` header row, sorted by
/// descending probability.
///
/// # Errors
/// Returns an error if the file cannot be created or written to.
pub fn write_results(
    path: impl AsRef<Path>,
    id_header: &str,
    results: HashMap<String, f32>,
) -> Result<(), Box<dyn Error>> {
    let mut rows: Vec<_> = results.into_iter().collect();
    rows.sort_by(|left, right| right.1.total_cmp(&left.1));

    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(&path)?;
    writer.write_record([id_header, "probability"])?;
    for (id, probability) in rows {
        writer.write_record([id, probability.to_string()])?;
    }
    writer.flush()?;
    println!("Results written to {}", fs::canonicalize(path)?.display());
    Ok(())
}

fn tsv_reader(path: impl AsRef<Path>) -> Result<csv::Reader<std::fs::File>, Box<dyn Error>> {
    Ok(ReaderBuilder::new().delimiter(b'\t').from_path(path)?)
}

/// Blocks the calling (synchronous) thread on `future` using a lazily-initialized Tokio runtime.
///
/// This CLI has no async call sites of its own — the only `async` code it reaches is
/// [`weight_effects::perform_effects_weighing_typed`]'s Unipept HTTP call for `taxonomic_analysis`
/// — so a single blocking entry point is simpler than threading `async`/`.await` through
/// `main()`, argument parsing, and every `read_*`/`write_results` helper. Mirrors the
/// `block_on_binding_future` helper `peptonizer_rust`'s PyO3 bindings use for the same reason.
fn block_on<F: Future>(future: F) -> F::Output {
    static TOKIO_RUNTIME: OnceLock<Runtime> = OnceLock::new();

    TOKIO_RUNTIME
        .get_or_init(|| {
            Builder::new_multi_thread()
                .enable_all()
                .build()
                .expect("Failed to initialize Tokio runtime")
        })
        .block_on(future)
}

/// Checks that every peptide referenced in `relationships` has a score, a count, and at least one
/// relationship ID.
///
/// # Errors
/// Returns an error describing the first missing score, missing count, empty ID list, or
/// completely empty `relationships` map that is found.
fn validate_inputs(
    relationships: &HashMap<String, Vec<usize>>,
    scores: &HashMap<String, f32>,
    counts: &HashMap<String, usize>,
) -> Result<(), Box<dyn Error>> {
    if relationships.is_empty() {
        return Err("The relationship TSV contains no rows".into());
    }

    for (peptide, relationships) in relationships {
        if !scores.contains_key(peptide) {
            return Err(format!("No score was provided for peptide {peptide}").into());
        }
        if !counts.contains_key(peptide) {
            return Err(format!("No count was provided for peptide {peptide}").into());
        }
        if relationships.is_empty() {
            return Err(format!("No IDs were provided for peptide {peptide}").into());
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::{Path, PathBuf};
    use std::time::{SystemTime, UNIX_EPOCH};

    /// A file in the OS temp directory that is deleted when dropped, so tests that write real
    /// files (the `read_*`/`write_results` functions all take `impl AsRef<Path>`) don't leak
    /// fixtures on disk, including when an assertion panics.
    struct TempFile(PathBuf);

    impl TempFile {
        fn with_content(name: &str, content: &str) -> Self {
            let nanos = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos();
            let path = std::env::temp_dir().join(format!("peptonizer_analysis_test_{name}_{nanos}"));
            fs::write(&path, content).expect("failed to write temp test fixture");
            TempFile(path)
        }

        fn empty(name: &str) -> Self {
            let nanos = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos();
            let path = std::env::temp_dir().join(format!("peptonizer_analysis_test_{name}_{nanos}"));
            TempFile(path)
        }

        fn path(&self) -> &Path {
            &self.0
        }
    }

    impl Drop for TempFile {
        fn drop(&mut self) {
            let _ = fs::remove_file(&self.0);
        }
    }

    // --- has_help_flag_in ---

    #[test]
    fn help_flag_detects_long_form() {
        let args = ["bin".to_string(), "--help".to_string()];
        assert!(has_help_flag_in(args.into_iter()));
    }

    #[test]
    fn help_flag_detects_short_form() {
        let args = ["bin".to_string(), "-h".to_string()];
        assert!(has_help_flag_in(args.into_iter()));
    }

    #[test]
    fn help_flag_absent_when_not_passed() {
        let args = ["bin".to_string(), "--peptide-taxa".to_string(), "f.tsv".to_string()];
        assert!(!has_help_flag_in(args.into_iter()));
    }

    // --- parse_arguments_from ---

    fn args(values: &[&str]) -> impl Iterator<Item = String> {
        values.iter().map(|value| value.to_string()).collect::<Vec<_>>().into_iter()
    }

    #[test]
    fn parse_arguments_happy_path() {
        let parsed = parse_arguments_from(
            args(&[
                "bin",
                "--peptide-taxa",
                "taxa.tsv",
                "--peptide-scores",
                "scores.tsv",
                "--peptide-counts",
                "counts.tsv",
                "--output",
                "out.tsv",
            ]),
            "--peptide-taxa",
            "default.tsv",
        )
        .unwrap();

        assert_eq!(parsed.relationships, "taxa.tsv");
        assert_eq!(parsed.scores, "scores.tsv");
        assert_eq!(parsed.counts, "counts.tsv");
        assert_eq!(parsed.output, "out.tsv");
    }

    #[test]
    fn parse_arguments_uses_default_output_when_omitted() {
        let parsed = parse_arguments_from(
            args(&[
                "bin",
                "--peptide-taxa",
                "taxa.tsv",
                "--peptide-scores",
                "scores.tsv",
                "--peptide-counts",
                "counts.tsv",
            ]),
            "--peptide-taxa",
            "default.tsv",
        )
        .unwrap();

        assert_eq!(parsed.output, "default.tsv");
    }

    #[test]
    fn parse_arguments_errors_on_missing_relationship_flag() {
        let result = parse_arguments_from(
            args(&["bin", "--peptide-scores", "s.tsv", "--peptide-counts", "c.tsv"]),
            "--peptide-taxa",
            "default.tsv",
        );
        assert!(result.is_err());
    }

    #[test]
    fn parse_arguments_errors_on_unknown_flag() {
        let result = parse_arguments_from(
            args(&["bin", "--totally-unknown", "value"]),
            "--peptide-taxa",
            "default.tsv",
        );
        assert!(result.is_err());
    }

    #[test]
    fn parse_arguments_errors_on_dangling_flag_with_no_value() {
        let result = parse_arguments_from(
            args(&["bin", "--peptide-taxa"]),
            "--peptide-taxa",
            "default.tsv",
        );
        assert!(result.is_err());
    }

    // --- read_relationships ---

    #[test]
    fn read_relationships_accumulates_multiple_ids_per_peptide() {
        let file = TempFile::with_content("relationships", "PEPTIDEA\t1\nPEPTIDEB\t2\nPEPTIDEA\t3\n");
        let relationships = read_relationships(file.path()).unwrap();
        assert_eq!(relationships["PEPTIDEA"], vec![1, 3]);
        assert_eq!(relationships["PEPTIDEB"], vec![2]);
    }

    #[test]
    fn read_relationships_skips_unparsable_header_row() {
        let file = TempFile::with_content("relationships_header", "peptide\tid\nPEPTIDEA\t1\n");
        let relationships = read_relationships(file.path()).unwrap();
        assert_eq!(relationships.len(), 1);
        assert_eq!(relationships["PEPTIDEA"], vec![1]);
    }

    #[test]
    fn read_relationships_errors_on_bad_id_past_first_row() {
        let file = TempFile::with_content("relationships_bad", "PEPTIDEA\t1\nPEPTIDEB\tnot_a_number\n");
        assert!(read_relationships(file.path()).is_err());
    }

    // --- read_relationships_with_string_ids ---

    #[test]
    fn read_relationships_with_string_ids_assigns_stable_ids_and_reverse_mapping() {
        let file = TempFile::with_content(
            "protein_relationships",
            "PEPTIDEA\tPROT1\nPEPTIDEB\tPROT2\nPEPTIDEC\tPROT1\n",
        );
        let (relationships, ids_by_index) =
            read_relationships_with_string_ids(file.path()).unwrap();

        let prot1_id = relationships["PEPTIDEA"][0];
        assert_eq!(relationships["PEPTIDEC"], vec![prot1_id]);
        assert_eq!(ids_by_index[&prot1_id], "PROT1");

        let prot2_id = relationships["PEPTIDEB"][0];
        assert_ne!(prot1_id, prot2_id);
        assert_eq!(ids_by_index[&prot2_id], "PROT2");
    }

    #[test]
    fn read_relationships_with_string_ids_accepts_non_numeric_ids() {
        let file = TempFile::with_content(
            "function_relationships",
            "PEPTIDEA\tGO:0006915\nPEPTIDEB\tEC:1.1.1.1\n",
        );
        let (relationships, ids_by_index) =
            read_relationships_with_string_ids(file.path()).unwrap();

        let go_id = relationships["PEPTIDEA"][0];
        assert_eq!(ids_by_index[&go_id], "GO:0006915");
    }

    // --- restore_original_ids ---

    #[test]
    fn restore_original_ids_maps_indices_back_to_ids() {
        let mut results = HashMap::new();
        results.insert("0".to_string(), 0.9f32);
        results.insert("1".to_string(), 0.1f32);

        let mut ids_by_index = HashMap::new();
        ids_by_index.insert(0usize, "PROT1".to_string());
        ids_by_index.insert(1usize, "PROT2".to_string());

        let restored = restore_original_ids(results, &ids_by_index).unwrap();
        assert_eq!(restored["PROT1"], 0.9f32);
        assert_eq!(restored["PROT2"], 0.1f32);
    }

    #[test]
    fn restore_original_ids_errors_on_unknown_index() {
        let mut results = HashMap::new();
        results.insert("42".to_string(), 0.5f32);

        let restored = restore_original_ids(results, &HashMap::new());
        assert!(restored.is_err());
    }

    // --- read_scores / read_counts ---

    #[test]
    fn read_scores_parses_values_and_skips_header() {
        let file = TempFile::with_content("scores", "peptide\tscore\nPEPTIDEA\t0.75\n");
        let scores = read_scores(file.path()).unwrap();
        assert_eq!(scores["PEPTIDEA"], 0.75f32);
    }

    #[test]
    fn read_counts_parses_values_and_skips_header() {
        let file = TempFile::with_content("counts", "peptide\tcount\nPEPTIDEA\t3\n");
        let counts = read_counts(file.path()).unwrap();
        assert_eq!(counts["PEPTIDEA"], 3);
    }

    // --- validate_inputs ---

    #[test]
    fn validate_inputs_accepts_consistent_data() {
        let mut relationships = HashMap::new();
        relationships.insert("PEPTIDEA".to_string(), vec![1]);
        let mut scores = HashMap::new();
        scores.insert("PEPTIDEA".to_string(), 0.5f32);
        let mut counts = HashMap::new();
        counts.insert("PEPTIDEA".to_string(), 1usize);

        assert!(validate_inputs(&relationships, &scores, &counts).is_ok());
    }

    #[test]
    fn validate_inputs_rejects_empty_relationships() {
        let result = validate_inputs(&HashMap::new(), &HashMap::new(), &HashMap::new());
        assert!(result.is_err());
    }

    #[test]
    fn validate_inputs_rejects_missing_score() {
        let mut relationships = HashMap::new();
        relationships.insert("PEPTIDEA".to_string(), vec![1]);
        let scores = HashMap::new();
        let mut counts = HashMap::new();
        counts.insert("PEPTIDEA".to_string(), 1usize);

        assert!(validate_inputs(&relationships, &scores, &counts).is_err());
    }

    #[test]
    fn validate_inputs_rejects_missing_count() {
        let mut relationships = HashMap::new();
        relationships.insert("PEPTIDEA".to_string(), vec![1]);
        let mut scores = HashMap::new();
        scores.insert("PEPTIDEA".to_string(), 0.5f32);
        let counts = HashMap::new();

        assert!(validate_inputs(&relationships, &scores, &counts).is_err());
    }

    #[test]
    fn validate_inputs_rejects_peptide_with_no_ids() {
        let mut relationships = HashMap::new();
        relationships.insert("PEPTIDEA".to_string(), vec![]);
        let mut scores = HashMap::new();
        scores.insert("PEPTIDEA".to_string(), 0.5f32);
        let mut counts = HashMap::new();
        counts.insert("PEPTIDEA".to_string(), 1usize);

        assert!(validate_inputs(&relationships, &scores, &counts).is_err());
    }

    // --- write_results ---

    #[test]
    fn write_results_sorts_descending_and_writes_header() {
        let mut results = HashMap::new();
        results.insert("low".to_string(), 0.1f32);
        results.insert("high".to_string(), 0.9f32);

        let out = TempFile::empty("write_results");
        write_results(out.path(), "protein_id", results).unwrap();

        let content = fs::read_to_string(out.path()).unwrap();
        let mut lines = content.lines();
        assert_eq!(lines.next().unwrap(), "protein_id\tprobability");
        assert!(lines.next().unwrap().starts_with("high\t"));
        assert!(lines.next().unwrap().starts_with("low\t"));
    }
}
