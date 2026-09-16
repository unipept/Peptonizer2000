extern crate serde_json;
extern crate serde;

mod utils;
mod http_client;
mod random;
pub mod weight_effects;
pub mod zero_lookahead_belief_propagation;
pub mod factor_graph;
mod fetch_unipept_taxa;
mod unipept_communicator;
pub mod effects_clustering;
pub mod analyse_grid_search;
#[cfg(not(target_arch = "wasm32"))]
mod input_parser;
#[cfg(not(target_arch = "wasm32"))]
mod clean_csv;

#[cfg(target_arch = "wasm32")]
pub use wasm::*;

#[cfg(not(target_arch = "wasm32"))]
pub use pyo3::*;

#[cfg(target_arch = "wasm32")]
mod wasm {
    use wasm_bindgen::prelude::*;
    use crate::fetch_unipept_taxa::fetch_peptides_and_filter_taxa;
    use crate::weight_effects::perform_effects_weighing;
    use crate::zero_lookahead_belief_propagation::run_belief_propagation;
    use crate::factor_graph::generate_graph;
    use crate::effects_clustering::cluster_effects;
    use crate::analyse_grid_search::compute_goodness;

    extern crate wasm_bindgen;
    extern crate web_sys;
    extern crate wasm_bindgen_futures;
    extern crate js_sys;
    extern crate console_error_panic_hook;

    /// Fetches taxa for peptides and filters them by rank and taxon query.
    ///
    /// # Arguments
    /// * `peptides` - JSON string of peptide sequences.
    /// * `rank` - Taxonomic rank used for filtering (e.g. "species").
    /// * `taxon_query` - JSON string of taxon IDs to filter against.
    ///
    /// # Returns
    /// JSON string mapping peptides to filtered taxon IDs.
    ///
    /// # Panics
    /// Panics if input JSON cannot be parsed or if result cannot be serialized.
    #[wasm_bindgen]
    pub async fn fetch_unipept_taxa_wasm(
        peptides: String,
        rank: String,
        taxon_query: String,
    ) -> Result<String, JsValue> {
        fetch_peptides_and_filter_taxa(peptides, rank, taxon_query)
            .await
            .map_err(|e| JsValue::from_str(&format!("fetch_unipept_taxa_wasm failed: {e}")))
    }

    /// Represents the main pipeline for weighting effects based on peptide evidence.
    ///
    /// # Arguments
    ///
    /// * `pep_effects` - JSON string mapping peptide sequences to lists of effect IDs.
    /// * `pep_scores` - JSON string mapping peptide sequences to their scores (float).
    /// * `pep_psm_counts` - JSON string mapping peptide sequences to their PSM counts (int).
    /// * `max_effects` - Maximum number of effects to include in output.
    /// * `effects_rank` - The effect rank to normalize effects to (e.g., "species"), or `None` to skip normalization.
    ///
    /// # Returns
    ///
    /// Tuple `(sequence_csv, effects_weights_csv)`:
    /// * `sequence_csv` - CSV string of peptide sequences and their weights.
    /// * `effects_weights_csv` - CSV string of effects weights and uniqueness.
    #[wasm_bindgen]
    pub async fn perform_effects_weighing_wasm(
        pep_effects: String,
        pep_scores: String,
        pep_psm_counts: String,
        max_effects: usize,
        effects_rank: Option<String>
    ) -> Result<Box<[JsValue]>, JsValue> {
        console_error_panic_hook::set_once(); // Enable panic logging
        let (sequence_csv, effects_weights_csv): (String, String) = perform_effects_weighing(pep_effects, pep_scores, pep_psm_counts, max_effects, effects_rank)
            .await
            .map_err(|e| JsValue::from_str(&format!("perform_effects_weighing_wasm failed: {e}")))?;

        Ok(Box::new([JsValue::from(sequence_csv), JsValue::from(effects_weights_csv)]))
    }

    /// Generates a GraphML representation of a factor graph from a CSV string of sequence scores.
    ///
    /// # Arguments
    /// * `sequence_scores_csv` - A string containing CSV data for sequence scores.
    ///
    /// # Returns
    /// Returns a `Result` containing a GraphML string representation of the factor graph.
    ///
    /// # Errors
    /// Returns an error if CSV parsing fails or if any error occurs during graph construction.
    #[wasm_bindgen]
    pub fn generate_pepgm_graph_wasm(sequence_scores_csv: String) -> Vec<u8> {
        let factor_graph_bytes = generate_graph(sequence_scores_csv).unwrap();

        factor_graph_bytes
    }

    /// Runs belief propagation on a factor graph provided as a GraphML string.
    ///
    /// This function constructs the factor graph, fills in factor tables and priors,
    /// splits the graph into connected components, and performs loopy belief propagation
    /// on each component. The result is returned as a CSV string.
    ///
    /// # Arguments
    ///
    /// * `graph` - GraphML representation of the factor graph.
    /// * `alpha` - Noisy-OR factor alpha parameter.
    /// * `beta` - Noisy-OR factor beta parameter.
    /// * `regularized` - Whether to regularize factor tables to penalize large numbers of parents.
    /// * `prior` - Prior belief for effect nodes.
    /// * `max_iter` - Maximum number of belief propagation iterations.
    /// * `tol` - Tolerance threshold for message convergence.
    ///
    /// # Returns
    ///
    /// CSV string with one row per node containing columns:
    /// `[node_name, posterior_probability_1, node_category]`
    #[wasm_bindgen]
    pub fn execute_pepgm_wasm(
        graphs: Vec<u8>,
        alpha: f32,
        beta: f32,
        regularized: bool,
        prior: f32,
        max_iter: Option<u32>,
        tol: Option<f32>
    ) -> String {
        // console_error_panic_hook::set_once(); // Enable panic logging

        run_belief_propagation(&graphs, alpha, beta, regularized, prior, max_iter, tol).unwrap()
    }

    /// Clusters effects based on peptidome similarity and returns a CSV.
    ///
    /// # Arguments
    /// * `sequence_scores_csv` - Sequence scores as CSV string.
    /// * `effects_weights_csv` - Effects weights as CSV string.
    /// * `similarity_threshold` - Threshold for clustering.
    ///
    /// # Returns
    /// CSV string with effects and their clusters.
    ///
    /// # Errors
    /// Returns an error if parsing, graph building, or clustering fails.
    #[wasm_bindgen]
    pub fn cluster_effects_wasm(
        sequence_scores_csv: String,
        effects_weights_csv: String,
        similarity_threshold: f32
    ) -> String {
        cluster_effects(sequence_scores_csv, effects_weights_csv, similarity_threshold).unwrap()
    }

    /// Computes a "goodness" score for clustering results by combining
    /// ranking similarity (via rank-biased overlap) and diversity (via entropy).
    /// 
    /// # Arguments
    /// * `effect_cluster_heads_csv` - CSV string file containing effect cluster heads.
    /// * `peptonizer_results` - JSON string containing effects scores produced by Peptonizer.
    /// 
    /// # Returns
    /// A `Result<f64, Box<dyn std::error::Error>>` containing the computed goodness score,
    /// or an error if parsing fails.
    /// 
    /// # Errors
    /// This function may return an error if the input CSV or JSON cannot be parsed.
    #[wasm_bindgen]
    pub fn compute_goodness_wasm(
        effect_cluster_heads_csv: String,
        peptonizer_results: String
    ) -> f64 {
        compute_goodness(&effect_cluster_heads_csv, &peptonizer_results).unwrap()
    }

}

#[allow(unsafe_op_in_unsafe_fn)]
#[cfg(not(target_arch = "wasm32"))]
mod pyo3 {
    use std::future::Future;
    use std::sync::OnceLock;
    use pyo3::prelude::*;
    use pyo3::types::PyBytes;
    use crate::fetch_unipept_taxa::fetch_peptides_and_filter_taxa;
    use crate::weight_effects::perform_effects_weighing;
    use crate::zero_lookahead_belief_propagation::run_belief_propagation;
    use crate::factor_graph::generate_graph;
    use crate::effects_clustering::cluster_effects;
    use crate::analyse_grid_search::compute_goodness;
    use crate::input_parser::{parse_input_peptides, parse_unique_peptides};
    use crate::clean_csv::clean_csv;
    use crate::unipept_communicator::get_names_for_taxa;

    extern crate console_error_panic_hook;

    fn block_on_binding_future<F>(future: F) -> F::Output
    where
        F: Future,
    {
        use tokio::runtime::{Builder, Runtime};

        static TOKIO_RUNTIME: OnceLock<Runtime> = OnceLock::new();

        TOKIO_RUNTIME
            .get_or_init(|| {
                Builder::new_multi_thread()
                    .enable_all()
                    .build()
                    .expect("Failed to initialize Tokio runtime for Python bindings")
            })
            .block_on(future)
    }

    /// Parses peptides from a TSV string and returns JSON representations
    /// of scores and counts.
    ///
    /// # Arguments
    /// * `tsv_content` - Input TSV string with peptide data.
    ///
    /// # Returns
    /// A tuple containing:
    /// * `String` - JSON of peptide → max score mapping.
    /// * `String` - JSON of peptide → occurrence count mapping.
    ///
    /// # Errors
    /// Returns an error if parsing fails or if JSON serialization fails.
    #[pyfunction]
    pub fn parse_input_peptides_py(tsv_content: String) -> (String, String) {
        parse_input_peptides(tsv_content).unwrap()
    }

    /// Extracts unique peptides from a TSV string and returns them as JSON.
    ///
    /// # Arguments
    /// * `tsv_content` - Input TSV string with peptide data.
    ///
    /// # Returns
    /// JSON string containing the list of unique peptides.
    ///
    /// # Errors
    /// Returns an error if parsing fails or if JSON serialization fails.
    #[pyfunction]
    pub fn parse_unique_peptides_py(tsv_content: String) -> String {
        parse_unique_peptides(tsv_content).unwrap()
    }

    /// Fetches effects for peptides and filters them by rank and effect query.
    ///
    /// # Arguments
    /// * `peptides` - JSON string of peptide sequences.
    /// * `rank` - Taxonomic rank used for filtering (e.g. "species").
    /// * `taxon_query` - JSON string of effect IDs to filter against.
    ///
    /// # Returns
    /// JSON string mapping peptides to filtered effect IDs.
    ///
    /// # Panics
    /// Panics if input JSON cannot be parsed or if result cannot be serialized.
    #[pyfunction]
    pub fn fetch_unipept_taxa_py(
        peptides: String,
        rank: String,
        taxon_query: String,
    ) -> String {
        block_on_binding_future(fetch_peptides_and_filter_taxa(
            peptides,
            rank,
            taxon_query,
        ))
        .unwrap()
    }

    /// Represents the main pipeline for weighting effects based on peptide evidence.
    ///
    /// # Arguments
    ///
    /// * `pep_effects` - JSON string mapping peptide sequences to lists of effect IDs.
    /// * `pep_scores` - JSON string mapping peptide sequences to their scores (float).
    /// * `pep_psm_counts` - JSON string mapping peptide sequences to their PSM counts (int).
    /// * `max_effects` - Maximum number of effects to include in output.
    /// * `effects_rank` - The effect rank to normalize effects to (e.g., "species").
    ///
    /// # Returns
    ///
    /// Tuple `(sequence_csv, effects_weights_csv)`:
    /// * `sequence_csv` - CSV string of peptide sequences and their weights.
    /// * `effects_weights_csv` - CSV string of effects weights and uniqueness.
    #[pyfunction]
    fn perform_effects_weighing_py(
        unipept_responses: String,
        pep_scores: String,
        pep_psm_counts: String,
        max_effects: usize,
        effects_rank: String
    ) -> (String, String) {
        block_on_binding_future(perform_effects_weighing(unipept_responses, pep_scores, pep_psm_counts, max_effects, Some(effects_rank))).unwrap()
    }

    /// Generates a GraphML representation of a factor graph from a CSV string of effect weights.
    ///
    /// # Arguments
    /// * `effects_weights_csv` - A string containing CSV data for effect weights.
    ///
    /// # Returns
    /// Returns a `Result` containing a GraphML string representation of the factor graph.
    ///
    /// # Errors
    /// Returns an error if CSV parsing fails or if any error occurs during graph construction.
    #[pyfunction]
    pub fn generate_pepgm_graph_py(py: Python<'_>, effects_weights_csv: String) -> Py<PyBytes> {
        let graph_bytes = generate_graph(effects_weights_csv).unwrap();

        PyBytes::new_bound(py, &graph_bytes).into()
    }

    /// Runs belief propagation on a factor graph provided as a GraphML string.
    ///
    /// This function constructs the factor graph, fills in factor tables and priors,
    /// splits the graph into connected components, and performs loopy belief propagation
    /// on each component. The result is returned as a CSV string.
    ///
    /// # Arguments
    ///
    /// * `graph` - GraphML representation of the factor graph.
    /// * `alpha` - Noisy-OR factor alpha parameter.
    /// * `beta` - Noisy-OR factor beta parameter.
    /// * `regularized` - Whether to regularize factor tables to penalize large numbers of parents.
    /// * `prior` - Prior belief for effect nodes.
    /// * `max_iter` - Maximum number of belief propagation iterations.
    /// * `tol` - Tolerance threshold for message convergence.
    ///
    /// # Returns
    ///
    /// CSV string with one row per node containing columns:
    /// `[node_name, posterior_probability_1, node_category]`
    #[pyfunction]
    #[pyo3(signature = (graph, alpha, beta, regularized, prior, max_iter=None, tol=None))]
    pub fn execute_pepgm_py(
        graph: Vec<u8>,
        alpha: f32,
        beta: f32,
        regularized: bool,
        prior: f32,
        max_iter: Option<u32>,
        tol: Option<f32>
    ) -> String {
        console_error_panic_hook::set_once(); // Enable panic logging

        run_belief_propagation(&graph, alpha, beta, regularized, prior, max_iter, tol).unwrap()
    }

    /// Clusters effects based on peptidome similarity and returns a CSV.
    ///
    /// # Arguments
    /// * `sequence_scores_csv` - CSV string containing peptide sequence scores.
    /// * `effects_weights_csv` - Effects weights as CSV string.
    /// * `similarity_threshold` - Threshold for clustering.
    ///
    /// # Returns
    /// CSV string with effects and their clusters.
    ///
    /// # Errors
    /// Returns an error if parsing, graph building, or clustering fails.
    #[pyfunction]
    pub fn cluster_effects_py(
        sequence_scores_csv: String,
        effects_weights_csv: String,
        similarity_threshold: f32
    ) -> String {
        cluster_effects(sequence_scores_csv, effects_weights_csv, similarity_threshold).unwrap()
    }

    /// Computes a "goodness" score for clustering results by combining
    /// ranking similarity (via rank-biased overlap) and diversity (via entropy).
    /// 
    /// # Arguments
    /// * `effect_cluster_heads_csv` - CSV string file containing effect cluster heads.
    /// * `peptonizer_results` - JSON string containing effects scores produced by Peptonizer.
    /// 
    /// # Returns
    /// A `Result<f64, Box<dyn std::error::Error>>` containing the computed goodness score,
    /// or an error if parsing fails.
    /// 
    /// # Errors
    /// This function may return an error if the input CSV or JSON cannot be parsed.
    #[allow(unsafe_op_in_unsafe_fn)]
    #[pyfunction]
    pub fn compute_goodness_py(
        effect_cluster_heads_csv: String,
        peptonizer_results: String
    ) -> f64 {
        compute_goodness(&effect_cluster_heads_csv, &peptonizer_results).unwrap()
    }

    /// Returns a mapping from effect ID to effect name for all effects provided.
    ///
    /// # Arguments
    /// * `target_effects` - A list of effect IDs for which all corresponding effect names should be retrieved.
    ///
    /// # Errors
    /// Returns an error if the Unipept API server responds with a non-success status code
    /// or if something goes wrong with the network or JSON parsing.
    ///
    /// # Returns
    /// A JSON string mapping effect IDs to their corresponding effect names.
    #[pyfunction]
    pub fn get_names_for_taxa_py(target_effects: Vec<usize>) -> String {
        let names = block_on_binding_future(get_names_for_taxa(&target_effects)).unwrap();
        serde_json::to_string(&names).unwrap()
    }

    /// Read a CSV-file that was produced by the PepGM algorithm and use it to
    /// produce a new CSV-file that only contains the effect-related information
    /// and scores. The string produced by this function can be written directly
    /// to a valid CSV-file and contains three columns: `effect_name`, `effect_id`,
    /// and `score`.
    ///
    /// # Arguments
    /// * `csv_content` - A CSV-file (as a string) that has been generated by running the PepGM algorithm.
    ///
    /// # Returns
    /// A `String` containing CSV rows with the columns: `effect_name,effect_id,score`.
    #[pyfunction]
    pub fn clean_csv_py(csv_content: String) -> String {
        block_on_binding_future(clean_csv(csv_content)).unwrap()
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[pymodule]
    fn peptonizer_rust(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(parse_input_peptides_py, m)?)?;
        m.add_function(wrap_pyfunction!(parse_unique_peptides_py, m)?)?;
        m.add_function(wrap_pyfunction!(fetch_unipept_taxa_py, m)?)?;
        m.add_function(wrap_pyfunction!(perform_effects_weighing_py, m)?)?;
        m.add_function(wrap_pyfunction!(generate_pepgm_graph_py, m)?)?;
        m.add_function(wrap_pyfunction!(execute_pepgm_py, m)?)?;
        m.add_function(wrap_pyfunction!(cluster_effects_py, m)?)?;
        m.add_function(wrap_pyfunction!(compute_goodness_py, m)?)?;
        m.add_function(wrap_pyfunction!(clean_csv_py, m)?)?;
        m.add_function(wrap_pyfunction!(get_names_for_taxa_py, m)?)?;
        Ok(())
    }
}
