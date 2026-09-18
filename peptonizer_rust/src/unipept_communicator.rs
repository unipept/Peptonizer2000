use std::collections::{HashMap, HashSet};
use std::cmp::min;
use crate::http_client::*;
use futures::future::join_all;
use serde::{Serialize, Deserialize};
use serde_json::{ Value };


/// Base URL for the UniPept API
const UNIPEPT_URL: &str = "https://api.unipept.ugent.be";
/// Endpoint for mapping peptides to filtered taxa
const UNIPEPT_PEPT2FILTERED_ENDPOINT: &str = "/api/v2/pept2taxa";
/// Endpoint for retrieving taxonomic lineages
const UNIPEPT_TAXONOMY_ENDPOINT: &str = "/api/v2/taxonomy";

/// Maximum number of peptides per request to the peptide-to-taxa endpoint
const UNIPEPT_PEPTIDES_BATCH_SIZE: usize = 2000;

/// Maximum number of taxa per request to the taxonomy endpoint
const TAXONOMY_ENDPOINT_BATCH_SIZE: usize = 100;


/// Standard NCBI taxonomy ranks for lineage retrieval
const NCBI_RANKS: &[&str] = &[
    "domain",
    "realm",
    "subkingdom",
    "superphylum",
    "phylum",
    "subphylum",
    "superclass",
    "class",
    "subclass",
    "superorder",
    "order",
    "suborder",
    "infraorder",
    "superfamily",
    "family",
    "subfamily",
    "tribe",
    "subtribe",
    "genus",
    "subgenus",
    "species_group",
    "species_subgroup",
    "species",
    "subspecies",
    "strain",
    "varietas",
    "forma"
];


/// Payload structure for requesting taxonomy lineages from UniPept
#[derive(Serialize, Deserialize, Debug)]
pub struct HTTPTaxonomyPayload {
    input: Vec<usize>,
    extra: bool
}

/// Payload structure for mapping peptides to taxa
#[derive(Serialize, Deserialize, Debug)]
pub struct HTTPPept2TaxaPayload { 
    input: Vec<String>,
    compact: bool,
    tryptic: bool
}

/// Response structure for peptide-to-taxa mapping
#[derive(Serialize, Deserialize, Debug)]
pub struct HTTPPept2TaxaResponse {
    peptide: String,
    taxa: Vec<usize>
}

/// Payload structure for retrieving descendants of taxa at specified ranks
#[derive(Serialize, Deserialize, Debug)]
pub struct HTTPTaxonomyDescendantsPayload {
    input: Vec<usize>,
    descendants: bool,
    descendants_ranks: Vec<String>
}

/// Response structure for retrieving descendants of a taxon
#[derive(Serialize, Deserialize, Debug)]
pub struct HTTPTaxonomyDescendantsResponse {
    taxon_id: usize,
    taxon_name: String,
    taxon_rank: String,
    descendants: Vec<usize>
}

/// Represents a response from the UniPept taxonomy API
#[cfg(not(target_arch = "wasm32"))]
#[derive(Serialize, Deserialize, Debug)]
struct TaxonomyResponse {
    taxon_id: usize,
    taxon_name: String,
}


/// Parses a JSON string returned by the UniPept API into a vector of key-value maps
///
/// # Arguments
/// * `http_response` - JSON string from UniPept API
///
/// # Returns
/// Vector of hash maps, where each key maps to an `Option<usize>` value. Only numeric or null values are retained.
#[allow(clippy::type_complexity)]
fn parse_response_json_string(http_response: &str) -> HttpResult<Vec<HashMap<String, Option<usize>>>> {
    let http_response_map: Vec<HashMap<String, Option<usize>>> = serde_json::from_str::<Vec<HashMap<String, Value>>>(http_response)?
            .into_iter()
            .map(|mut obj: HashMap<String, Value>| {
                // Remove the key-value pair where the value is a string
                obj.retain(|_, v| v.is_null() || v.is_number());

                // Convert the remaining keys and values to `HashMap<String, Option<usize>>`
                obj.into_iter()
                    .map(|(key, value)| {
                        let value = if value.is_number() {
                            Some(value.as_i64().unwrap() as usize)
                        } else {
                            None
                        };
                        (key, value)
                    })
                    .collect()
            })
            .collect();
    
    Ok(http_response_map)
}

/// Retrieves the unique lineage taxa IDs at a specified taxonomic rank.
///
/// This function queries the UniPept taxonomy API for the given `target_taxa` and extracts
/// the taxon IDs at the specified `taxa_rank`. To minimize API requests, it uses a cache
/// (`lineage_cache`) to store previously fetched lineages. 
///
/// # Arguments
///
/// * `target_taxa` - A reference to a vector of taxon IDs for which the lineage is requested.
/// * `taxa_rank` - The target taxonomic rank (e.g., "species", "genus") at which the unique lineage is extracted.
/// * `lineage_cache` - A mutable reference to a hash map that stores previously fetched lineages.
///
/// # Returns
///
/// A vector of unique taxon IDs corresponding to the specified `taxa_rank`.
///
/// # Panics
///
/// The function will panic if:
/// - The `taxa_rank` does not exist in the predefined `NCBI_RANKS`.
pub async fn get_unique_lineage_at_specified_rank_async(target_taxa: &[usize], taxa_rank: &str, lineage_cache: &mut HashMap<usize, Vec<Option<usize>>>) -> HttpResult<Vec<usize>> {

    let url: String = [UNIPEPT_URL, UNIPEPT_TAXONOMY_ENDPOINT].concat();

    // Remove duplicates from input
    let target_taxa: HashSet<usize> = target_taxa.iter().cloned().collect();

    // Prepare a list of taxa that are not yet in the cache
    let taxa_to_request: Vec<usize> = target_taxa.iter().filter(| tax_id | ! lineage_cache.contains_key(tax_id)).cloned().collect();

    let http_client = &create_http_client();
    let mut pending_requests = Vec::new();

    for i in (0..taxa_to_request.len()).step_by(TAXONOMY_ENDPOINT_BATCH_SIZE) {

        let batch_size: usize = std::cmp::min(TAXONOMY_ENDPOINT_BATCH_SIZE, taxa_to_request.len() - i);
        let batch: Vec<usize> = taxa_to_request[i..(i + batch_size)].to_vec();
        let payload = HTTPTaxonomyPayload { input: batch, extra: true };
        let payload_json = serde_json::to_string(&payload)?;

        pending_requests.push(http_client.perform_post_request(url.clone(), payload_json));
    }

    let responses = join_all(pending_requests).await;

    // Parse responses one-by-one after all requests have completed.
    for (batch_idx, http_response) in responses.into_iter().enumerate() {
        let http_response: String = http_response
            .map_err(|e| format!("Failed to retrieve taxonomy data for batch {}. Error message: {}", batch_idx, e))?;
        let http_response = parse_response_json_string(&http_response)?;

        for lineage_json in http_response {
            let lineage_json: HashMap<String, Option<usize>> = lineage_json;
            let lineage: Vec<Option<usize>> = NCBI_RANKS.iter()
                    .filter_map(|key| lineage_json.get(&format!("{key}_id")).cloned())
                    .collect();
            let taxon_id: usize = lineage_json.get("taxon_id").ok_or("Taxon ID not in lineage")?.ok_or("Taxon ID is None")?;
            lineage_cache.insert(taxon_id, lineage);
        }

    }

    let rank_idx = NCBI_RANKS.iter().position(|&ncbi_rank| ncbi_rank == taxa_rank).ok_or("Taxa rank not found in NCBI ranks")?;
    let lineage: HashSet<usize> = target_taxa.iter()
                                            .filter_map(|taxon| lineage_cache.get(taxon).and_then(|lineage| lineage[rank_idx]))
                                            .collect();
    let lineage: Vec<usize> = lineage.into_iter().collect();

    Ok(lineage)
}

/// Queries Unipept and returns all the taxa that are associated with the given list of peptides.
/// 
/// For each peptide in the input, an entry in the output map is created, which points to the
/// taxon IDs associated with this peptide.
/// 
/// # Arguments
/// 
/// * `peptides` - A list of peptide sequences for which all associated taxa should be queried.
/// 
/// # Errors
/// 
/// Returns an error if the Unipept API server responds with an error or if a network issue occurs.
/// 
/// # Returns
/// 
/// A map from each peptide in the input list to its associated taxa IDs.
pub async fn get_taxa_for_peptides_async(peptides: Vec<String>) -> HttpResult<HashMap<String, Vec<usize>>> {
    
    let url = [UNIPEPT_URL, UNIPEPT_PEPT2FILTERED_ENDPOINT].concat();

    let mut output = HashMap::new();

    let http_client = &create_http_client();
    let mut pending_requests = Vec::new();

    for i in (0..peptides.len()).step_by(UNIPEPT_PEPTIDES_BATCH_SIZE) {

        let end_batch = min(i+UNIPEPT_PEPTIDES_BATCH_SIZE, peptides.len());
        let batch = peptides[i..end_batch].to_vec();
        let payload = HTTPPept2TaxaPayload { input: batch, compact: true, tryptic: true };
        let payload_json = serde_json::to_string(&payload)?;

        pending_requests.push(http_client.perform_post_request(url.clone(), payload_json));
    }

    let responses = join_all(pending_requests).await;

    // Parse responses one-by-one after all requests have completed.
    for (batch_idx, http_response) in responses.into_iter().enumerate() {
        let http_response: String = http_response
            .map_err(|e| format!("Failed to retrieve taxa data for batch {}. Error message: {}", batch_idx, e))?;

        let http_response = serde_json::from_str::<Vec<HTTPPept2TaxaResponse>>(&http_response)?;

        for peptide_data in &http_response {
            let original_taxa: Vec<usize> = peptide_data.taxa.clone();
            output.insert(peptide_data.peptide.clone(), original_taxa);
        }
    }

    Ok(output)
}

/// Returns a list of all taxon IDs that are descendants of the given taxa in `target_taxa`.
///
/// # Arguments
///
/// * `target_taxa` - A list of taxon IDs for which all descendants at a specific NCBI rank (and lower) should be retrieved.
/// * `descendants_rank` - The maximum rank that each of the descendants should have in the NCBI taxonomy.
///   All descendants that are defined at this rank or deeper are reported.
///
/// # Errors
///
/// Returns an error if the Unipept API server responds with an error, or if something goes wrong
/// with the network.
///
/// # Returns
///
/// A list of taxon IDs that meet the given rank criteria.
pub async fn get_descendants_for_taxa_async(target_taxa: Vec<usize>, descendant_rank: String) -> HttpResult<HashSet<usize>> {
    
    let url = [UNIPEPT_URL, UNIPEPT_TAXONOMY_ENDPOINT].concat();
    let mut all_descendants = HashSet::new();

    // We need to get all children at the requested level, AND at lower levels. That's what we are using the ranks array for.
    let rank_idx = NCBI_RANKS.iter().position(|&ncbi_rank| ncbi_rank == descendant_rank).ok_or("descendants rank not found in NCBI ranks")?;
    let descentants_ranks: Vec<String> = NCBI_RANKS[rank_idx..].iter().map(|&s| s.to_string()).collect();

    let http_client = &create_http_client();
    let mut pending_requests = Vec::new();

    for i in (0..target_taxa.len()).step_by(TAXONOMY_ENDPOINT_BATCH_SIZE) {

        let end_batch = min(i+TAXONOMY_ENDPOINT_BATCH_SIZE, target_taxa.len());
        let batch = target_taxa[i..end_batch].to_vec();
        let payload = HTTPTaxonomyDescendantsPayload { input: batch, descendants: true, descendants_ranks: descentants_ranks.clone() };
        let payload_json = serde_json::to_string(&payload)?;

        pending_requests.push(http_client.perform_post_request(url.clone(), payload_json));
    }

    let responses = join_all(pending_requests).await;

    // Parse responses one-by-one after all requests have completed.
    for (batch_idx, http_response) in responses.into_iter().enumerate() {
        let http_response = http_response
            .map_err(|e| format!("Failed to retrieve taxonomy data for batch {}. Error message: {}", batch_idx, e))?;
        let http_response = serde_json::from_str::<Vec<HTTPTaxonomyDescendantsResponse>>(&http_response)?;
        
        for response in http_response {
            all_descendants.extend(response.descendants);
        }
    }

    Ok(all_descendants)
}

/// Returns a mapping from taxon ID to taxon name for all taxa provided.
///
/// # Arguments
/// * `target_taxa` - A list of taxon IDs for which all corresponding taxon names should be retrieved.
///
/// # Errors
/// Returns an error if the Unipept API server responds with a non-success status code
/// or if something goes wrong with the network or JSON parsing.
///
/// # Returns
/// A `HashMap<usize, String>` mapping taxon IDs to their corresponding taxon names.
#[cfg(not(target_arch = "wasm32"))]
pub async fn get_names_for_taxa(target_taxa: &[usize]) -> HttpResult<HashMap<usize, String>> {
    let url = format!("{UNIPEPT_URL}{UNIPEPT_TAXONOMY_ENDPOINT}");
    let mut output: HashMap<usize, String> = HashMap::new();

    let http_client = &create_http_client();
    let mut pending_requests = Vec::new();

    for i in (0..target_taxa.len()).step_by(TAXONOMY_ENDPOINT_BATCH_SIZE) {
        let batch: Vec<usize> = target_taxa[i..std::cmp::min(i + TAXONOMY_ENDPOINT_BATCH_SIZE, target_taxa.len())]
            .to_vec();

        let payload = serde_json::json!({
            "input": batch
        });
        let payload_json = serde_json::to_string(&payload)?;

        pending_requests.push(http_client.perform_post_request(url.clone(), payload_json));
    }

    let responses = join_all(pending_requests).await;

    // Parse responses one-by-one after all requests have completed.
    for http_response in responses {
        let http_response = http_response.map_err(|e| format!("Communication error: {e}"))?;

        let http_response = serde_json::from_str::<Vec<TaxonomyResponse>>(&http_response)?;
        
        for response in http_response {
            output.insert(response.taxon_id, response.taxon_name);
        }
    }

    Ok(output)
}



#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_parse_response_json_string() {
        let http_response = r#"
        [
            {"taxon_id": 1, "species_id": 9606, "genus_id": null, "name": "Homo sapiens"},
            {"taxon_id": 2, "species_id": null, "genus_id": 9605, "name": "Pan troglodytes"}
        ]
        "#;

        let parsed = parse_response_json_string(http_response);
        assert!(parsed.is_ok());
        let parsed = parsed.unwrap();

        assert_eq!(parsed.len(), 2);
        assert_eq!(parsed[0].get("taxon_id"), Some(&Some(1)));
        assert_eq!(parsed[0].get("species_id"), Some(&Some(9606)));
        assert_eq!(parsed[0].get("genus_id"), Some(&None));
        assert!(!parsed[0].contains_key("name"));
        assert_eq!(parsed[1].get("taxon_id"), Some(&Some(2)));
        assert_eq!(parsed[1].get("genus_id"), Some(&Some(9605)));
    }

    #[tokio::test]
    async fn test_get_unique_lineage_at_specified_rank_with_cache() {
        let mut lineage_cache: HashMap<usize, Vec<Option<usize>>> = HashMap::new();
        lineage_cache.insert(1, vec![Some(1), Some(10), Some(100)]);
        lineage_cache.insert(2, vec![Some(2), Some(20), Some(200)]);

        let target_taxa = vec![1, 2];
        let rank = NCBI_RANKS[1];
        let lineage = get_unique_lineage_at_specified_rank_async(&target_taxa, rank, &mut lineage_cache).await;
        assert!(lineage.is_ok());
        let lineage = lineage.unwrap();

        assert!(lineage.contains(&10));
        assert!(lineage.contains(&20));
        assert_eq!(lineage.len(), 2);
    }

    #[tokio::test]
    async fn test_get_descendants_for_taxa_structure() {
        let descendants = get_descendants_for_taxa_async(vec![200, 701], "species".to_string()).await;
        assert!(descendants.is_ok());
        let descendants = descendants.unwrap();

        assert!(descendants.len() == 4);
    }

    #[tokio::test]
    async fn test_get_names_for_taxa_structure() {
        let taxa = vec![1, 2];
        let result = get_names_for_taxa(&taxa).await;
        assert!(result.is_ok());
        
        let names = result.unwrap();
        assert_eq!(names.get(&1).unwrap(), "root");
        assert_eq!(names.get(&2).unwrap(), "Bacteria");
    }

    #[tokio::test]
    async fn test_get_taxa_for_peptides_structure() {
        let peptides = vec!["AAEEAAAA".to_string(), "AAAAEEA".to_string()];
        let result = get_taxa_for_peptides_async(peptides).await;
        assert!(result.is_ok());
        let result = result.unwrap();

        // Check structure (both peptides present, each mapped to at least one taxon)
        assert!(!result.get("AAEEAAAA").unwrap().is_empty());
        assert!(!result.get("AAAAEEA").unwrap().is_empty());
    }
}
