use std::collections::HashMap;
use nori::zero_lookahead_bp_from_graph_bytes;


/// Runs belief propagation on a serialized factor graph (as produced by `generate_graph`).
///
/// This function constructs the factor graph, fills in factor tables and priors,
/// splits the graph into connected components, and performs loopy belief propagation
/// on each component. The result is returned as a JSON string.
///
/// # Arguments
///
/// * `graph` - Serialized factor graph bytes.
/// * `alpha` - Noisy-OR factor alpha parameter.
/// * `beta` - Noisy-OR factor beta parameter.
/// * `regularized` - Whether to regularize factor tables to penalize large numbers of parents.
/// * `prior` - Prior belief for effect nodes.
/// * `max_iter` - Maximum number of belief propagation iterations.
/// * `tol` - Tolerance threshold for message convergence.
///
/// # Returns
///
/// JSON string mapping each node name to its posterior probability.
pub fn run_belief_propagation(
    graph_bytes: &[u8],
    alpha: f32,
    beta: f32,
    regularized: bool,
    prior: f32,
    max_iter: Option<u32>,
    tol: Option<f32>
) -> Result<String, Box<dyn std::error::Error>> {
    let results = zero_lookahead_bp_from_graph_bytes(graph_bytes, alpha, beta, regularized, prior, max_iter, tol).unwrap();

    let effect_score_dict: HashMap<String, f32> = results
        .into_iter()
        .filter_map(|(key, values)| {
            values.get(1).map(|&v| (key, v))
        })
        .collect();

    Ok(serde_json::to_string(&effect_score_dict)?)
}

