use std::collections::{HashSet, HashMap};
use crate::effects_clustering::{Effect, parse_effect_csv};


/// Computes a "goodness" score for clustering results by combining
/// ranking similarity (via rank-biased overlap) and diversity (via entropy).
/// 
/// # Arguments
/// * `effect_cluster_heads_csv` - CSV string containing effect cluster heads.
/// * `peptonizer_results` - JSON string containing effects scores produced by Peptonizer.
/// 
/// # Returns
/// A `Result<f64, Box<dyn std::error::Error>>` containing the computed goodness score,
/// or an error if parsing fails.
/// 
/// # Errors
/// This function may return an error if the input CSV or JSON cannot be parsed.
pub fn compute_goodness(
    effect_cluster_heads_csv: &str,
    peptonizer_results: &str
) -> Result<f64, Box<dyn std::error::Error>> {

    let taxid_weights: Vec<Effect> = parse_effect_csv(effect_cluster_heads_csv)?;
    let effect: Vec<usize> = taxid_weights.iter().map(|t| t.effect).collect();

    let effects_scores: HashMap<String, f64> = serde_json::from_str(peptonizer_results)?;
    let mut effects_scores: Vec<(&String, &f64)> = effects_scores.iter().collect();

    effects_scores.sort_by(|a, b| b.1.partial_cmp(a.1).expect("Partial compare returned None")); // ascending order

    let sorted_ids: Vec<usize> = effects_scores.iter()
        .map(|(k, _)| (*k).clone().parse::<usize>())
        .collect::<Result<Vec<_>, _>>()?;
    let sorted_scores: Vec<f64> = effects_scores.iter().map(|(_, v)| **v).collect::<Vec<_>>();

    let entropy = entropy(&sorted_scores);
    let rbo = rbo(&effect, &sorted_ids);

    Ok(rbo / entropy.powi(2))
}


/// Computes the Rank-Biased Overlap (RBO) between two ranked lists of effect IDs.
/// RBO measures the agreement between two ranked lists, emphasizing higher ranks.
/// 
/// # Arguments
/// * `list1` - First ranked list of effect IDs.
/// * `list2` - Second ranked list of effect IDs.
/// 
/// # Returns
/// A value between 0.0 and 1.0 representing the similarity of the two ranked lists.
fn rbo(list1: &[usize], list2: &[usize]) -> f64 {
    let k = list1.len().min(list2.len());
    let mut sum = 0.0;

    let mut seen1 = HashSet::new();
    let mut seen2 = HashSet::new();

    for d in 1..=k {
        seen1.insert(list1[d - 1]);
        seen2.insert(list2[d - 1]);

        let overlap = seen1.intersection(&seen2).count() as f64;
        let agreement = overlap / d as f64;

        sum += agreement;
    }

    sum / k as f64
}


/// Calculates the Shannon entropy of a set of values.
/// Entropy measures the diversity or unpredictability of a distribution.
/// 
/// # Arguments
/// * `values` - A slice of floating-point values representing weights or probabilities.
/// 
/// # Returns
/// The entropy as a floating-point value. Returns 0.0 if all values sum to zero.
fn entropy(values: &[f64]) -> f64 {
    let sum: f64 = values.iter().sum();

    if sum == 0.0 {
        return 0.0;
    }

    values.iter()
        .map(|&v| {
            let p = v / sum;
            if p > 0.0 {
                -p * p.log2()
            } else {
                0.0
            }
        })
        .sum()
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::effects_clustering::generate_effects_cluster_csv;

    #[test]
    fn test_entropy_uniform_distribution() {
        let values = vec![1.0, 1.0, 1.0, 1.0];
        let result = entropy(&values);
        // Uniform distribution of 4 elements has entropy = log2(4) = 2.0
        assert!((result - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_entropy_all_zero() {
        let values = vec![0.0, 0.0, 0.0];
        let result = entropy(&values);
        assert_eq!(result, 0.0);
    }

    #[test]
    fn test_rbo_perfect_match() {
        let list1 = vec![1, 2, 3, 4];
        let list2 = vec![1, 2, 3, 4];
        let result = rbo(&list1, &list2);
        assert!((result - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_rbo_no_overlap() {
        let list1 = vec![1, 2, 3];
        let list2 = vec![4, 5, 6];
        let result = rbo(&list1, &list2);
        assert_eq!(result, 0.0);
    }

    #[test]
    fn test_compute_goodness_valid_inputs() {
        // Prepare a small CSV with two effects
        let effects = vec![
            Effect {
                id: 0,
                effect: 1,
                scaled_weight: 0.5,
                unique: true,
            },
            Effect {
                id: 1,
                effect: 2,
                scaled_weight: 0.8,
                unique: false,
            },
        ];
        let csv = generate_effects_cluster_csv(effects).unwrap();

        // JSON scores
        let json_scores = serde_json::json!({
            "1": 0.9,
            "2": 0.8
        }).to_string();

        let result = compute_goodness(&csv, &json_scores);
        assert!(result.is_ok());
        let score = result.unwrap();
        assert!(score.is_finite());
        assert!(score > 0.0);
    }
}