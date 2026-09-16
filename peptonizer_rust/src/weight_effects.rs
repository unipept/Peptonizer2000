use std::collections::{HashMap, HashSet};
use crate::utils::log;
use crate::random::select_random_samples_with_weights;
use csv::Writer;
use crate::unipept_communicator::get_unique_lineage_at_specified_rank_async;

/// Represents the main pipeline for weighting effects based on peptide evidence.
///
/// # Arguments
///
/// * `pep_effects` - JSON string mapping peptide peptides to lists of effect IDs.
/// * `pep_scores` - JSON string mapping peptide peptides to their scores (float).
/// * `pep_psm_counts` - JSON string mapping peptide peptides to their PSM counts (int).
/// * `max_effects` - Maximum number of effects to include in output.
/// * `effects_rank` - NCBI rank at which the Peptonizer analysis should be performed. Should be a rank that is supported by Unipept.
///
/// # Returns
///
/// Tuple `(sequence_csv, effects_weights_csv)`:
/// * `sequence_csv` - CSV string of peptide peptides and their weights.
/// * `effects_weights_csv` - CSV string of effects weights and uniqueness.
pub async fn perform_effects_weighing(
    pep_effects: String,
    pep_scores: String,
    pep_psm_counts: String,
    max_effects: usize,
    effects_rank: Option<String>
) -> Result<(String, String), Box<dyn std::error::Error>> {
    log("Parsing Unipept responses from disk...");
    let peptide_effects: HashMap<String, Vec<usize>> = serde_json::from_str(&pep_effects)?;
    let peptide_scores_map: HashMap<String, f32> = serde_json::from_str(&pep_scores)?;
    let peptide_counts_map: HashMap<String, usize> = serde_json::from_str(&pep_psm_counts)?;

    perform_effects_weighing_typed(peptide_effects, peptide_scores_map, peptide_counts_map, max_effects, effects_rank).await
}

/// Same pipeline as [`perform_effects_weighing`], but for native Rust callers (e.g. `peptonizer_analysis`)
/// that already have the peptide relationships, scores, and counts as typed maps and don't need the
/// JSON serialize/deserialize round-trip required by the WASM/PyO3 bindings.
///
/// # Arguments
///
/// * `peptide_effects` - Mapping of peptide sequences to lists of effect IDs.
/// * `peptide_scores_map` - Mapping of peptide sequences to their scores (float).
/// * `peptide_counts_map` - Mapping of peptide sequences to their PSM counts (int).
/// * `max_effects` - Maximum number of effects to include in output.
/// * `effects_rank` - NCBI rank at which the Peptonizer analysis should be performed. Should be a rank that is supported by Unipept.
///
/// # Returns
///
/// Tuple `(sequence_csv, effects_weights_csv)`:
/// * `sequence_csv` - CSV string of peptide peptides and their weights.
/// * `effects_weights_csv` - CSV string of effects weights and uniqueness.
pub async fn perform_effects_weighing_typed(
    peptide_effects: HashMap<String, Vec<usize>>,
    peptide_scores_map: HashMap<String, f32>,
    peptide_counts_map: HashMap<String, usize>,
    max_effects: usize,
    effects_rank: Option<String>
) -> Result<(String, String), Box<dyn std::error::Error>> {
    let peptides: Vec<String> = peptide_effects.keys().map(|seq| seq.to_owned()).collect();

    let mut effects: Vec<Vec<usize>> = peptide_effects.into_values().collect();

    if let Some(ref rank) = effects_rank {
        log("Started mapping all effect ids to the specified rank...");
        normalize_unipept_responses(&mut effects, rank).await?;
    } else {
        log("Skipping rank normalization because no rank was provided...");
    }
    let chosen_effects: HashSet<usize> = weighted_random_sample(&effects, 10000)?;

    log(&format!("Using {} peptides as input...", chosen_effects.len()));

    log("Normalizing peptides and converting to vectors...");

    let peptides: Vec<String> = chosen_effects.iter().map(|idx| peptides[*idx].to_owned()).collect();
    let effects: Vec<Vec<usize>> = chosen_effects.iter().map(|idx| effects[*idx].to_owned()).collect();

    // Only keep the randomly selected samples.
    let mut peptide_scores: Vec<f32> = vec![0.0; peptides.len()];
    for i in 0..peptides.len() {
        peptide_scores[i] = peptide_scores_map[&peptides[i]];
    }

    let mut peptide_counts: Vec<usize> = vec![0; peptides.len()];
    for i in 0..peptides.len() {
        peptide_counts[i] = peptide_counts_map[&peptides[i]];
    }

    /* Score the degeneracy of a effects, i.e.,
       how conserved a peptide sequence is between effects.
       map all taxids in the list in the effects column back to their taxid at species level (or the rank specified by the user)
       Right now, Effect is simply a copy of effects. This step still needs to be optimized.
    */

    // Divide the number of PSMs of a peptide by the number of effects the peptide is associated with, exponentiated by 3
    log("Started dividing the number of PSMS of a peptide by the number the peptide is associated with...");
    let peptide_weights: Vec<f32> = peptide_counts.iter()
                                          .zip(effects.iter().map(|connected_effects| connected_effects.len().pow(3)))
                                          .map(|(&count, len_cube)| count as f32 / len_cube as f32)
                                          .collect();

    let unique_effects: HashSet<usize> = effects.iter()
                                                   .filter(|tax| tax.len() == 1)
                                                   .map(|tax| tax[0])
                                                   .collect();

    // Sum up the weights of a effect and sort by weight
    log("Started summing the weights of a effect and sorting them by weight...");
    let peptide_log_weights: Vec<f32> = peptide_weights.iter().map(|w| (w + 1.0).log10()).collect();

    //  Since large proteomes tend to have more detectable peptides,
    // we adjust the weight by dividing by the size of the proteome i.e.,
    // the number of proteins that are associated with a effect
    let peptide_scaled_weight = peptide_log_weights.clone();

    let mut effect_weights: HashMap<usize, f32> = HashMap::new();
    for (ids, weight) in effects.clone().into_iter().zip(peptide_scaled_weight.clone().into_iter()) {
        for id in ids {
            *effect_weights.entry(id).or_insert(0.0) += weight;
        }
    }
    let mut sorted_effect_weights: Vec<(usize, f32)> = effect_weights.into_iter().collect();
    sorted_effect_weights.sort_by(|a, b| b.1.partial_cmp(&a.1).expect("Partial compare returned None"));
    let (effects_sorted, effect_weights_sorted): (Vec<usize>, Vec<f32>) = sorted_effect_weights.into_iter().unzip();

    let unique_effects_vec: Vec<bool> = effects_sorted.iter().map(|id| unique_effects.contains(id)).collect();

    let peptides_csv = if effect_weights_sorted.len() < 50 {
        generate_peptides_csv(None, false, peptides, peptide_scores, peptide_counts, effects, peptide_weights, peptide_log_weights)?
    } else {
        let mut effects_to_include: HashSet<usize> = effects_sorted.iter().take(max_effects).cloned().collect();
        effects_to_include.extend(unique_effects);

        generate_peptides_csv(Some(effects_to_include), true, peptides, peptide_scores, peptide_counts, effects, peptide_weights, peptide_log_weights)?
    };

    let effects_weights_csv = generate_effects_weights_csv(effects_sorted, effect_weights_sorted, unique_effects_vec)?;

    Ok((peptides_csv, effects_weights_csv))
}

/// Generates a CSV for peptides with associated effect weights and scores.
///
/// # Arguments
///
/// * `effects_to_include` - Optional set of effects IDs to filter the output.
/// * `filter_effects` - Whether to filter peptides based on `effects_to_include`.
/// * `peptides` - List of peptide sequences.
/// * `scores` - List of peptide scores corresponding to `peptides`.
/// * `psms` - List of PSM counts corresponding to `peptides`.
/// * `effect` - List of lists of higher effects IDs for each peptide.
/// * `weights` - Computed weights for each peptide.
/// * `log_weights` - Log-transformed weights for each peptide.
///
/// # Returns
///
/// CSV string containing one row per peptide-effect pair with columns:
/// "id", "sequence", "score", "psms", "effect", "weight", "log_weight".
#[allow(clippy::too_many_arguments)]
fn generate_peptides_csv(effects_to_include: Option<HashSet<usize>>, filter_effects: bool, peptides: Vec<String>, scores: Vec<f32>, psms: Vec<usize>, effect: Vec<Vec<usize>>, weights: Vec<f32>, log_weights: Vec<f32>) -> Result<String, Box<dyn std::error::Error>> {

    let mut wtr = Writer::from_writer(vec![]);

    let _ = wtr.write_record(["id", "sequence", "score", "psms", "effect", "weight", "log_weight"]);

    let mut id = 0;
    for i in 0..peptides.len() {
        for effect in &effect[i] {
            if (! filter_effects) || effects_to_include.as_ref().ok_or("No effects to include passed while filter_effects enabled")?.contains(effect) {
                wtr.write_record(&[
                    id.to_string(),
                    peptides[i].clone(), 
                    scores[i].to_string(), 
                    psms[i].to_string(), 
                    effect.to_string(), 
                    weights[i].to_string(), 
                    log_weights[i].to_string()
                ])?;
                id += 1;
            }
        }
    }

    let csv: String = String::from_utf8(wtr.into_inner()?)?;

    Ok(csv)
}

/// Generates a CSV of effects weights.
///
/// # Arguments
///
/// * `effect` - List of effect IDs.
/// * `higher_taxid_weights` - List of computed weights corresponding to `effect`.
/// * `higher_taxid_unique` - List indicating whether each effect is uniquely associated with a peptide.
///
/// # Returns
///
/// CSV string with columns: "id", "effect", "scaled_weight", "unique".
fn generate_effects_weights_csv(effect: Vec<usize>, higher_taxid_weights: Vec<f32>, higher_taxid_unique: Vec<bool>) -> Result<String, Box<dyn std::error::Error>> {
    let mut wtr = Writer::from_writer(vec![]);

    let _ = wtr.write_record(["id", "effect", "scaled_weight", "unique"]);

    for i in 0..effect.len() {
        wtr.write_record(&[
            i.to_string(),
            effect[i].to_string(),
            higher_taxid_weights[i].to_string(),
            higher_taxid_unique[i].to_string()
        ])?;
    }

    let csv: String = String::from_utf8(wtr.into_inner()?)?;

    Ok(csv)
    
}

/// Maps effects lists onto the effect rank specified by the user.
///
/// # Arguments
///
/// * `effects` - Mutable reference to a vector of vectors of effect IDs.
/// * `effects_rank` - The desired effect rank to normalize to (e.g., "species").
pub(crate) async fn normalize_unipept_responses(effects: &mut [Vec<usize>], effects_rank: &str) -> Result<(), Box<dyn std::error::Error>> {

    // TODO: should we first do get_lineages_for_effects to limit Unipept calls (see python)?
    let mut lineage_cache: HashMap<usize, Vec<Option<usize>>> = HashMap::new();

    // Map all effects onto the rank specified by the user
    for effect in effects {
        *effect = get_unique_lineage_at_specified_rank_async(effect, effects_rank, &mut lineage_cache)
            .await
            .map_err(|e| -> Box<dyn std::error::Error> { e })?;
    }

    Ok(())
}

/// Selects `n` random indices from `effects` vectors, weighted by inverse degeneracy.
///
/// # Arguments
///
/// * `effects` - Vector of vectors of effect IDs for each peptide.
/// * `n` - Number of random samples to select.
///
/// # Returns
///
/// A `HashSet` of selected indices, chosen with probability proportional to
/// `1 / number_of_effects_per_peptide`.
fn weighted_random_sample(effects: &[Vec<usize>], n: usize) -> Result<HashSet<usize>, Box<dyn std::error::Error>> {
    
    // Calculate normalized weights based on the length of the effects array
    let weights: Vec<f64> = effects.iter().map(|effect| if effect.is_empty() { 0.0 } else { 1.0 / effect.len() as f64 }).collect();
    let total_weight: f64 = weights.iter().sum();
    let normalized_weights: Vec<f64> = weights.iter().map(|w| w / total_weight).collect();

    let samples: HashSet<usize> = select_random_samples_with_weights(normalized_weights, n)?;

    Ok(samples)
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[tokio::test]
    async fn test_perform_effects_weighing_basic() {
        let pep_effects_json = r#"{"PEP1":[3000],"PEP2":[3500]}"#.to_string();
        let pep_scores_json = r#"{"PEP1":0.8,"PEP2":0.5}"#.to_string();
        let pep_psm_counts_json = r#"{"PEP1":4,"PEP2":2}"#.to_string();
        let max_effects = 10;
        let effects_rank = None;

        let csvs = perform_effects_weighing(
            pep_effects_json,
            pep_scores_json,
            pep_psm_counts_json,
            max_effects,
            effects_rank
        )
        .await;
        assert!(csvs.is_ok());
        let (seq_csv, effects_csv) = csvs.unwrap();

        assert!(seq_csv.contains("sequence"));
        assert!(seq_csv.contains("PEP1"));
        assert!(effects_csv.contains("effect"));
    }

    #[test]
    fn test_generate_sequence_csv_basic() {
        let peptides = vec!["PEP1".to_string(), "PEP2".to_string()];
        let scores = vec![0.8, 0.5];
        let psms = vec![4, 2];
        let effect = vec![vec![3000], vec![3001]];
        let weights = vec![0.5, 0.2];
        let log_weights = vec![0.18, 0.079];

        let csv = generate_peptides_csv(None, false, peptides, scores, psms, effect, weights, log_weights);
        assert!(csv.is_ok());
        let csv = csv.unwrap();
        assert!(csv.contains("sequence"));
        assert!(csv.contains("PEP1"));
        assert!(csv.contains("1"));
    }

    #[test]
    fn test_generate_sequence_csv_with_filter() {
        let peptides = vec!["PEP".to_string()];
        let scores = vec![0.9];
        let psms = vec![2];
        let effect = vec![vec![10,11,12]];
        let weights = vec![0.2];
        let log_weights = vec![0.042];
        let filter_effects: HashSet<usize> = vec![11,12].into_iter().collect();

        let csv = generate_peptides_csv(Some(filter_effects), true, peptides, scores, psms, effect, weights, log_weights);
        assert!(csv.is_ok());
        let csv = csv.unwrap();
        assert!(csv.contains("12"));
        assert!(!csv.contains("10")); 
    }

    #[test]
    fn test_generate_effects_weights_csv_basic() {
        let effect = vec![1,2];
        let weights = vec![0.5, 0.8];
        let unique_flags = vec![true, false];

        let csv = generate_effects_weights_csv(effect, weights, unique_flags);
        assert!(csv.is_ok());
        let csv = csv.unwrap();
        assert!(csv.contains("effect"));
        assert!(csv.contains("0")); 
        assert!(csv.contains("true"));
    }

    #[tokio::test]
    async fn test_normalize_unipept_responses_basic() {
        let mut effects = vec![vec![3000]];
        let _ = normalize_unipept_responses(&mut effects, "species").await;
        assert!(effects.iter().all(|v| !v.is_empty()));
    }

    #[test]
    fn test_weighted_random_sample_basic() {
        let effects = vec![vec![1], vec![2,3], vec![4]];
        let n = 2;
        let samples = weighted_random_sample(&effects, n);
        assert!(samples.is_ok());
        let samples = samples.unwrap();

        assert_eq!(samples.len(), n);
        assert!(samples.iter().all(|&idx| idx < effects.len()));
    }
}
