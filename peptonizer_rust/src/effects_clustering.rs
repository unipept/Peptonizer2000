use std::collections::{HashMap, HashSet};
use serde::{Serialize, Deserialize};
use csv::{ReaderBuilder, WriterBuilder};
use crate::factor_graph::{parse_effect_weights_csv, EffectWeight};


/// Represents a effect unit with attributes used for clustering.
#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct Effect {
    /// Unique identifier of the effect.
    pub id: usize,
    /// Identifier of the effect it belongs to.
    pub effect: usize,
    /// Weight of the effect, scaled for comparison.
    pub scaled_weight: f32,
    /// Whether the effect is unique in the dataset.
    pub unique: bool
}


/// Parses a CSV string into a list of `Effect`.
///
/// # Arguments
/// * `effects_weights_csv` - CSV content as string.
///
/// # Returns
/// Vector of `Effect`.
///
/// # Errors
/// Returns an error if CSV parsing fails.
pub fn parse_effect_csv(effects_weights_csv: &str) -> Result<Vec<Effect>, Box<dyn std::error::Error>> {
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .from_reader(effects_weights_csv.as_bytes());
    
    let mut effects_weights = Vec::new();
    for record in rdr.deserialize() {
        let row: Effect = record?;
        effects_weights.push(row);
    }

    Ok(effects_weights)
}


/// Serializes a list of `Effect` into CSV format.
///
/// # Arguments
/// * `effects` - List of effects to serialize.
///
/// # Returns
/// CSV as a string.
///
/// # Errors
/// Returns an error if serialization fails.
pub fn generate_effects_cluster_csv(effects: Vec<Effect>) -> Result<String, Box<dyn std::error::Error>> {
    let mut wtr = WriterBuilder::new().from_writer(vec![]);

    for effect in &effects {
        wtr.serialize(effect)?;
    }

    let data = String::from_utf8(wtr.into_inner()?)?;

    Ok(data)
}


/// Cluster effects based on peptidome similarity and returns a CSV containing the effects that are cluster heads.
///
/// # Arguments
/// * `graph_xml` - GraphML as string.
/// * `effects_weights_csv` - Effects weights as CSV string.
/// * `similarity_threshold` - Threshold for clustering.
///
/// # Returns
/// CSV string with effects and their clusters.
///
/// # Errors
/// Returns an error if parsing, graph building, or clustering fails.
pub fn cluster_effects(sequence_scores_csv: String, effects_weights_csv: String, similarity_threshold: f32) -> Result<String, Box<dyn std::error::Error>> {

    let sequence_scores = parse_effect_weights_csv(sequence_scores_csv)?;
    let effects_weights = parse_effect_csv(&effects_weights_csv)?;

    let peptidome_dict = get_peptides_per_effect(&sequence_scores)?;
    let (similarities, effect_index) = compute_detected_peptidome_similarity(peptidome_dict);

    let effects_weights_sorted: Vec<Effect> = effects_weights
        .into_iter()
        .filter(|tw| effect_index.contains_key(&tw.effect))
        .collect();
    let mut effects_sorted: Vec<usize> = effects_weights_sorted.iter().map(|tw| tw.effect).collect();
    let mut cluster_heads: Vec<usize> = Vec::new();

    while ! effects_sorted.is_empty() {
        let effect1 = effects_sorted[0];
        let mut cluster_list: Vec<usize> = Vec::new();
        cluster_heads.push(effect1);

        for &effect2 in &effects_sorted {
            if similarities[effect_index[&effect2]][effect_index[&effect1]] > similarity_threshold {
                cluster_list.push(effect2);
            }
        }
        effects_sorted.retain(|&effect| ! cluster_list.contains(&effect));

    }

    let effect_weights_sorted: Vec<Effect> = effects_weights_sorted.into_iter()
        .filter(|tw| cluster_heads.contains(&tw.effect))
        .collect();

    let effect_cluster_heads_csv = generate_effects_cluster_csv(effect_weights_sorted)?;
    Ok(effect_cluster_heads_csv)
}


/// Builds a dictionary of peptides per effect.
///
/// # Arguments
/// * `graph` - Factor graph reference.
///
/// # Returns
/// Map from effect ID to peptide set.
///
/// # Errors
/// Returns an error if node parsing fails.
fn get_peptides_per_effect(effect_weights: &Vec<EffectWeight>) -> Result<HashMap<usize, HashSet<usize>>, Box<dyn std::error::Error>> {
    let mut peptidome_dict = HashMap::new();

    // maps unique sequence string -> generated sequence_id
    let mut sequence_to_id: HashMap<String, usize> = HashMap::new();
    let mut next_sequence_id = 0usize;
    for tw in effect_weights {
        // get or create sequence_id
        let sequence_id = match sequence_to_id.get(&tw.sequence) {
            Some(id) => *id,
            None => {
                let id = next_sequence_id;
                next_sequence_id += 1;
                sequence_to_id.insert(tw.sequence.clone(), id);
                id
            }
        };
        // insert sequence_id into effect set
        peptidome_dict
            .entry(tw.effect)
            .or_insert_with(HashSet::new)
            .insert(sequence_id);
    }

    Ok(peptidome_dict)
}


/// Computes similarity matrix and effect index map.
///
/// # Arguments
/// * `peptidome_dict` - Map of effect to peptides.
///
/// # Returns
/// Tuple of (similarity matrix, effect index map).
fn compute_detected_peptidome_similarity(peptidome_dict: HashMap<usize, HashSet<usize>>) -> (Vec<Vec<f32>>, HashMap<usize, usize>) {
    let mut sim_matrix = Vec::new();
    let mut effect_index: HashMap<usize, usize> = HashMap::new();

    let peptidome_keys = peptidome_dict.keys();
    for (index, effect1) in peptidome_keys.clone().enumerate() {
        effect_index.insert(*effect1, index);
        let set1 = &peptidome_dict[effect1];
        let mut sim_row = Vec::new();
        for effect2 in peptidome_keys.clone() {
            let set2 = &peptidome_dict[effect2];
            let shared = set1.intersection(set2).count();
            let sim: f32 = if set2.is_empty() {
                0.0
            } else {
                shared as f32 / set2.len() as f32
            };

            sim_row.push(sim);
        }
        sim_matrix.push(sim_row);
    }

    (sim_matrix, effect_index)
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::{HashMap, HashSet};

    #[test]
    fn test_vec_to_string_and_string_to_vec() {
        let values = vec![1, 2, 3];
        let serialized = serde_json::to_string(&values).unwrap();
        assert_eq!(serialized, "[1,2,3]");

        let deserialized: Vec<usize> = serde_json::from_str(&serialized).unwrap();
        assert_eq!(deserialized, values);
    }

    #[test]
    fn test_parse_effect_csv_and_generate_effects_cluster_csv() {
        let csv_data = "\
id,effect,scaled_weight,unique
1,10,0.5,true
2,20,0.8,false
";

        let effects = parse_effect_csv(csv_data).unwrap();
        assert_eq!(effects.len(), 2);
        assert_eq!(effects[0].id, 1);
        assert_eq!(effects[1].effect, 20);

        let out_csv = generate_effects_cluster_csv(effects).unwrap();
        assert!(out_csv.contains("id,effect,scaled_weight,unique"));
        assert!(out_csv.contains("10"));
    }

    #[test]
    fn test_compute_detected_peptidome_similarity() {
        let mut peptidome_dict: HashMap<usize, HashSet<usize>> = HashMap::new();
        peptidome_dict.insert(1, HashSet::from([1, 2]));
        peptidome_dict.insert(2, HashSet::from([2, 3]));

        let (sim_matrix, effect_index) = compute_detected_peptidome_similarity(peptidome_dict);
        println!("{:?}", sim_matrix);

        assert_eq!(sim_matrix.len(), 2);
        assert_eq!(sim_matrix[0].len(), 2);
        assert!(effect_index.contains_key(&1));
        assert!(effect_index.contains_key(&2));
        assert!(sim_matrix[0][1] > 0.0);
    }

    #[test]
    fn test_generate_effects_cluster_csv_roundtrip() {
        let effects = vec![
            Effect { id: 1, effect: 10, scaled_weight: 0.5, unique: true },
            Effect { id: 2, effect: 20, scaled_weight: 0.8, unique: false },
        ];

        let csv_string = generate_effects_cluster_csv(effects.clone()).unwrap();
        let parsed = parse_effect_csv(&csv_string).unwrap();

        assert_eq!(parsed.len(), 2);
        assert_eq!(parsed[0].effect, 10);
        assert_eq!(parsed[1].effect, 20);
    }
}
