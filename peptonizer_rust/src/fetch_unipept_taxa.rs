use crate::unipept_communicator::{get_taxa_for_peptides_async, get_descendants_for_taxa_async};
use crate::http_client::HttpResult;
use std::collections::{HashMap, HashSet};


/// Fetches taxa for peptides and keeps only the taxa that are descendants of the taxon query.
///
/// The returned taxa are **not** normalized to `rank`: `rank` is only used to select which
/// descendants of `taxon_query` are kept. Normalization to `rank` is done by
/// [`crate::weight_effects::perform_effects_weighing`] when it receives an `effects_rank`, or
/// must be done by the caller before weighing (as `peptonizer_ts` expects).
///
/// # Arguments
/// * `peptides` - JSON string of peptide sequences.
/// * `rank` - Highest rank a descendant of `taxon_query` may have to be kept (e.g. "species").
/// * `taxon_query` - JSON string of taxon IDs to filter against.
///
/// # Returns
/// JSON string mapping peptides to filtered taxon IDs.
///
/// # Panics
/// Panics if input JSON cannot be parsed or if result cannot be serialized.
pub async fn fetch_peptides_and_filter_taxa(
    peptides: String,
    rank: String,
    taxon_query: String
) -> HttpResult<String> {
    // Parse arguments
    let peptides: Vec<String> = serde_json::from_str(&peptides)?;
    let taxon_query_ids: Vec<usize> = serde_json::from_str(&taxon_query)?;
    
    // First we retrieve all taxa associated with the given peptids
    let mut peptides_taxa: HashMap<String, Vec<usize>> = get_taxa_for_peptides_async(peptides).await?;

    // Then, we make sure to filter the taxa and only keep those that are associated 
    // to the taxa of interest indicated by the user. Retrieve all (in)direct children
    // of the filter taxa provided by the user
    let taxa_filter: HashSet<usize> = get_descendants_for_taxa_async(taxon_query_ids, rank.clone()).await?;

    // Compute the intersection of the taxa that should be retained and the original list of taxa
    for taxa_list in peptides_taxa.values_mut() {
        taxa_list.retain(|taxon| taxa_filter.contains(taxon));
    }

    Ok(serde_json::to_string(&peptides_taxa)?)
}


#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;

    #[tokio::test]
    async fn test_fetch_with_known_peptide_and_species() {
        let peptides = serde_json::to_string(&vec!["TATAAAA".to_string()]).unwrap();

        let taxon_query = serde_json::to_string(&vec![2]).unwrap();

        let result = fetch_peptides_and_filter_taxa(peptides, "species".to_string(), taxon_query).await;
        assert!(result.is_ok());
        let result = result.unwrap();

        let parsed: Value = serde_json::from_str(&result).unwrap();
        assert!(parsed.is_object());

        assert!(parsed.get("TATAAAA").is_some());
    }

    #[tokio::test]
    async fn test_empty_peptides_and_taxa() {
        let peptides = "[]".to_string();
        let taxon_query = "[]".to_string();

        let result = fetch_peptides_and_filter_taxa(peptides, "species".to_string(), taxon_query).await;
        assert!(result.is_ok());
        let result = result.unwrap();

        let parsed: Value = serde_json::from_str(&result).unwrap();
        assert!(parsed.is_object());

        assert_eq!(parsed.as_object().unwrap().len(), 0);
    }
}
