use std::collections::{HashMap, HashSet};
use crate::utils::log;
use crate::random::select_random_samples_with_weights;
use csv::Writer;
use crate::unipept_communicator::get_unique_lineage_at_specified_rank_async;

/// Represents the main pipeline for weighting taxa based on peptide evidence.
///
/// # Arguments
///
/// * `pep_taxa` - JSON string mapping peptide sequences to lists of taxon IDs.
/// * `pep_scores` - JSON string mapping peptide sequences to their scores (float).
/// * `pep_psm_counts` - JSON string mapping peptide sequences to their PSM counts (int).
/// * `max_taxa` - Maximum number of taxa to include in output.
/// * `taxa_rank` - Optional NCBI rank at which taxa should be normalized.
///   If `None`, input taxa are assumed to already be normalized.
///
/// # Returns
///
/// Tuple `(sequence_csv, taxa_weights_csv)`:
/// * `sequence_csv` - CSV string of peptide sequences and their weights.
/// * `taxa_weights_csv` - CSV string of taxa weights and uniqueness.
pub async fn perform_taxa_weighing(
    pep_taxa: String,
    pep_scores: String,
    pep_psm_counts: String,
    max_taxa: usize,
    taxa_rank: Option<String>
) -> Result<(String, String), Box<dyn std::error::Error>> {
    log("Parsing Unipept responses from disk...");
    let pep_taxa: HashMap<String, Vec<usize>> = serde_json::from_str(&pep_taxa)?;

    let sequences: Vec<String> = pep_taxa.keys().map(|seq| seq.to_owned()).collect();

    let mut taxa: Vec<Vec<usize>> = pep_taxa.into_values().collect();
    
    if let Some(taxa_rank) = taxa_rank {
        log("Started mapping all taxon ids to the specified rank...");
        normalize_unipept_responses(&mut taxa, &taxa_rank).await?;
    } else {
        log("Skipping taxon normalization because taxa_rank is None...");
    }

    let chosen_idx: HashSet<usize> = weighted_random_sample(&taxa, 10000)?;

    log(&format!("Using {} sequences as input...", chosen_idx.len()));

    log("Normalizing peptides and converting to vectors...");

    let sequences: Vec<String> = chosen_idx.iter().map(|idx| sequences[*idx].to_owned()).collect();
    let taxa: Vec<Vec<usize>> = chosen_idx.iter().map(|idx| taxa[*idx].to_owned()).collect();

    // Parse scores from JSON string to hashmap, only keep the randomly selected samples.
    let pep_scores_map: HashMap<String, f32> = serde_json::from_str(&pep_scores)?;
    let mut pep_scores: Vec<f32> = vec![0.0; sequences.len()];
    for i in 0..sequences.len() {
        pep_scores[i] = pep_scores_map[&sequences[i]];
    }

    // parse counts from JSON string to hashmap, only keep the randomly selected samples.
    let pep_psm_counts_map: HashMap<String, usize> = serde_json::from_str(&pep_psm_counts)?;
    let mut pep_psm_counts: Vec<usize> = vec![0; sequences.len()];
    for i in 0..sequences.len() {
        pep_psm_counts[i] = pep_psm_counts_map[&sequences[i]];
    }

    /* Score the degeneracy of a taxa, i.e.,
       how conserved a peptide sequence is between taxa.
       map all taxids in the list in the taxa column back to their taxid at species level (or the rank specified by the user)
       Right now, HigherTaxa is simply a copy of taxa. This step still needs to be optimized.
       Move taxa to highertaxa because taxa is not used anymore.
    */
    let higher_taxa: Vec<Vec<usize>> = taxa; 

    // Divide the number of PSMs of a peptide by the number of taxa the peptide is associated with, exponentiated by 3
    log("Started dividing the number of PSMS of a peptide by the number the peptide is associated with...");
    let weights: Vec<f32> = pep_psm_counts.iter()
                                          .zip(higher_taxa.iter().map(|taxa| taxa.len().pow(3)))
                                          .map(|(&count, len_cube)| count as f32 / len_cube as f32)
                                          .collect();

    let unique_psm_taxa: HashSet<usize> = higher_taxa.iter()
                                                   .filter(|tax| tax.len() == 1)
                                                   .map(|tax| tax[0])
                                                   .collect();

    // Sum up the weights of a taxon and sort by weight
    log("Started summing the weights of a taxon and sorting them by weight...");
    let log_weights: Vec<f32> = weights.iter().map(|w| (w + 1.0).log10()).collect();

    //  Since large proteomes tend to have more detectable peptides,
    // we adjust the weight by dividing by the size of the proteome i.e.,
    // the number of proteins that are associated with a taxon
    let scaled_weight = log_weights.clone();

    let mut tax_id_weights: HashMap<usize, f32> = HashMap::new();
    for (ids, weight) in higher_taxa.clone().into_iter().zip(scaled_weight.clone()) {
        for id in ids {
            *tax_id_weights.entry(id).or_insert(0.0) += weight;
        }
    }
    let mut sorted_tax_id_weights: Vec<(usize, f32)> = tax_id_weights.into_iter().collect();
    sorted_tax_id_weights.sort_by(|a, b| b.1.partial_cmp(&a.1).expect("Partial compare returned None"));
    let (tax_ids, tax_id_weights): (Vec<usize>, Vec<f32>) = sorted_tax_id_weights.into_iter().unzip();

    // Retrieves the specified taxonomic rank taxid in the lineage of each of the species-level taxids returned by
    // Unipept for both the UnipeptFrame and the TaxIdWeightFrame
    let higher_unique_psm_taxids = unique_psm_taxa;

    // Group the duplicate entries of higher up taxa and sum their weights
    let higher_taxid_weights = tax_id_weights;

    let higher_taxid_unique: Vec<bool> = tax_ids.iter().map(|id| higher_unique_psm_taxids.contains(id)).collect();

    // TODO: Why hardcoded < 50
    let sequence_csv = if higher_taxid_weights.len() < 50 {
        generate_sequence_csv(None, false, sequences, pep_scores, pep_psm_counts, higher_taxa, weights, log_weights)?
    } else {
        let mut taxa_to_include: HashSet<usize> = tax_ids.iter().take(max_taxa).cloned().collect();
        taxa_to_include.extend(higher_unique_psm_taxids);

        generate_sequence_csv(Some(taxa_to_include), true, sequences, pep_scores, pep_psm_counts, higher_taxa, weights, log_weights)?
    };

    let taxa_weights_csv = generate_taxa_weights_csv(tax_ids, higher_taxid_weights, higher_taxid_unique)?;

    Ok((sequence_csv, taxa_weights_csv))
}

/// Generates a CSV for sequences with associated taxonomic weights and scores.
///
/// # Arguments
///
/// * `taxa_to_include` - Optional set of taxa IDs to filter the output.
/// * `filter_taxa` - Whether to filter sequences based on `taxa_to_include`.
/// * `sequences` - List of peptide sequences.
/// * `scores` - List of peptide scores corresponding to `sequences`.
/// * `psms` - List of PSM counts corresponding to `sequences`.
/// * `higher_taxa` - List of lists of higher taxa IDs for each sequence.
/// * `weights` - Computed weights for each sequence.
/// * `log_weights` - Log-transformed weights for each sequence.
///
/// # Returns
///
/// CSV string containing one row per peptide-taxon pair with columns:
/// "id", "sequence", "score", "psms", "higher_taxa", "weight", "log_weight".
#[allow(clippy::too_many_arguments)]
fn generate_sequence_csv(taxa_to_include: Option<HashSet<usize>>, filter_taxa: bool, sequences: Vec<String>, scores: Vec<f32>, psms: Vec<usize>, higher_taxa: Vec<Vec<usize>>, weights: Vec<f32>, log_weights: Vec<f32>) -> Result<String, Box<dyn std::error::Error>> {

    let mut wtr = Writer::from_writer(vec![]);

    let _ = wtr.write_record(["id", "sequence", "score", "psms", "higher_taxa", "weight", "log_weight"]);

    let mut id = 0;
    for i in 0..sequences.len() {
        for taxon in &higher_taxa[i] {
            if (! filter_taxa) || taxa_to_include.as_ref().ok_or("No taxa to include passed while filter_taxa enabled")?.contains(taxon) {
                wtr.write_record(&[
                    id.to_string(),
                    sequences[i].clone(), 
                    scores[i].to_string(), 
                    psms[i].to_string(), 
                    taxon.to_string(), 
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

/// Generates a CSV of taxa weights.
///
/// # Arguments
///
/// * `higher_taxa` - List of taxon IDs.
/// * `higher_taxid_weights` - List of computed weights corresponding to `higher_taxa`.
/// * `higher_taxid_unique` - List indicating whether each taxon is uniquely associated with a peptide.
///
/// # Returns
///
/// CSV string with columns: "id", "higher_taxa", "scaled_weight", "unique".
fn generate_taxa_weights_csv(higher_taxa: Vec<usize>, higher_taxid_weights: Vec<f32>, higher_taxid_unique: Vec<bool>) -> Result<String, Box<dyn std::error::Error>> {
    let mut wtr = Writer::from_writer(vec![]);

    let _ = wtr.write_record(["id", "higher_taxa", "scaled_weight", "unique"]);

    for i in 0..higher_taxa.len() {
        wtr.write_record(&[
            i.to_string(),
            higher_taxa[i].to_string(),
            higher_taxid_weights[i].to_string(),
            higher_taxid_unique[i].to_string()
        ])?;
    }

    let csv: String = String::from_utf8(wtr.into_inner()?)?;

    Ok(csv)
    
}

/// Maps taxa lists onto the taxonomic rank specified by the user.
///
/// # Arguments
///
/// * `taxa` - Mutable reference to a vector of vectors of taxon IDs.
/// * `taxa_rank` - The desired taxonomic rank to normalize to (e.g., "species").
pub(crate) async fn normalize_unipept_responses(taxa: &mut [Vec<usize>], taxa_rank: &str) -> Result<(), Box<dyn std::error::Error>> {
    
    // TODO: should we first do get_lineages_for_taxa to limit Unipept calls (see python)?
    let mut lineage_cache: HashMap<usize, Vec<Option<usize>>> = HashMap::new();

    // Map all taxa onto the rank specified by the user
    for taxon in taxa {
        *taxon = get_unique_lineage_at_specified_rank_async(taxon, taxa_rank, &mut lineage_cache)
            .await
            .map_err(|e| -> Box<dyn std::error::Error> { e })?;
    }

    Ok(())
}

/// Selects `n` random indices from `taxa` vectors, weighted by inverse degeneracy.
///
/// # Arguments
///
/// * `taxa` - Vector of vectors of taxon IDs for each peptide.
/// * `n` - Number of random samples to select.
///
/// # Returns
///
/// A `HashSet` of selected indices, chosen with probability proportional to
/// `1 / number_of_taxa_per_peptide`.
fn weighted_random_sample(taxa: &[Vec<usize>], n: usize) -> Result<HashSet<usize>, Box<dyn std::error::Error>> {
    
    // Calculate normalized weights based on the length of the taxa array
    let weights: Vec<f64> = taxa.iter().map(|taxon| if taxon.is_empty() { 0.0 } else { 1.0 / taxon.len() as f64 }).collect();
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
    async fn test_perform_taxa_weighing_basic() {
        let pep_taxa_json = r#"{"PEP1":[3000],"PEP2":[3500]}"#.to_string();
        let pep_scores_json = r#"{"PEP1":0.8,"PEP2":0.5}"#.to_string();
        let pep_psm_counts_json = r#"{"PEP1":4,"PEP2":2}"#.to_string();
        let max_taxa = 10;
        let taxa_rank = None;

        let csvs = perform_taxa_weighing(
            pep_taxa_json,
            pep_scores_json,
            pep_psm_counts_json,
            max_taxa,
            taxa_rank
        )
        .await;
        assert!(csvs.is_ok());
        let (seq_csv, taxa_csv) = csvs.unwrap();

        assert!(seq_csv.contains("sequence"));
        assert!(seq_csv.contains("PEP1"));
        assert!(taxa_csv.contains("higher_taxa"));
    }

    #[test]
    fn test_generate_sequence_csv_basic() {
        let sequences = vec!["PEP1".to_string(), "PEP2".to_string()];
        let scores = vec![0.8, 0.5];
        let psms = vec![4, 2];
        let higher_taxa = vec![vec![3000], vec![3001]];
        let weights = vec![0.5, 0.2];
        let log_weights = vec![0.18, 0.079];

        let csv = generate_sequence_csv(None, false, sequences, scores, psms, higher_taxa, weights, log_weights);
        assert!(csv.is_ok());
        let csv = csv.unwrap();
        assert!(csv.contains("sequence"));
        assert!(csv.contains("PEP1"));
        assert!(csv.contains("1"));
    }

    #[test]
    fn test_generate_sequence_csv_with_filter() {
        let sequences = vec!["PEP".to_string()];
        let scores = vec![0.9];
        let psms = vec![2];
        let higher_taxa = vec![vec![10,11,12]];
        let weights = vec![0.2];
        let log_weights = vec![0.042];
        let filter_taxa: HashSet<usize> = vec![11,12].into_iter().collect();

        let csv = generate_sequence_csv(Some(filter_taxa), true, sequences, scores, psms, higher_taxa, weights, log_weights);
        assert!(csv.is_ok());
        let csv = csv.unwrap();
        assert!(csv.contains("12"));
        assert!(!csv.contains("10")); 
    }

    #[test]
    fn test_generate_taxa_weights_csv_basic() {
        let higher_taxa = vec![1,2];
        let weights = vec![0.5, 0.8];
        let unique_flags = vec![true, false];

        let csv = generate_taxa_weights_csv(higher_taxa, weights, unique_flags);
        assert!(csv.is_ok());
        let csv = csv.unwrap();
        assert!(csv.contains("higher_taxa"));
        assert!(csv.contains("0")); 
        assert!(csv.contains("true"));
    }

    #[tokio::test]
    async fn test_normalize_unipept_responses_basic() {
        let mut taxa = vec![vec![3000]];
        let _ = normalize_unipept_responses(&mut taxa, "species").await;
        assert!(taxa.iter().all(|v| !v.is_empty()));
    }

    #[test]
    fn test_weighted_random_sample_basic() {
        let taxa = vec![vec![1], vec![2,3], vec![4]];
        let n = 2;
        let samples = weighted_random_sample(&taxa, n);
        assert!(samples.is_ok());
        let samples = samples.unwrap();

        assert_eq!(samples.len(), n);
        assert!(samples.iter().all(|&idx| idx < taxa.len()));
    }
}
