
use std::collections::HashSet;


/// Selects `n` unique random sample indices from a list of weights.
///
/// The selection is based on weighted probabilities, where higher weights increase
/// the likelihood of being chosen. Sampling is without replacement (indices are unique).
///
/// # Arguments
/// * `weights` - Vector of non-negative weights for each element.
/// * `n` - Number of unique indices to select.
///
/// # Returns
/// A `HashSet<usize>` containing the selected indices.
///
/// # Errors
/// This function will panic if:
/// * `weights` is empty.
/// * `n` is larger than `weights.len()`.
// #[cfg(not(target_arch = "wasm32"))]
pub fn select_random_samples_with_weights(
    weights: Vec<f64>,
    n: usize,
) -> Result<HashSet<usize>, Box<dyn std::error::Error>> {
    use rand::prelude::*;
    let mut rng = thread_rng();

    let mut keys: Vec<(f64, usize)> = weights.iter().enumerate()
        .filter(|&(_, &w)| w > 0.0)
        .map(|(i, &w)| {
            let u: f64 = rng.gen_range(0.0..1.0);
            let key = u.powf(1.0/w);
            (key, i)
        }).collect();
    
    keys.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    let selected_idx: HashSet<usize> = keys.iter()
        .take(n.min(keys.len()))
        .map(|&(_, i)| i)
        .collect();

    Ok(selected_idx)
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn test_basic_sampling() {
        let weights = vec![1.0, 2.0, 3.0];
        let samples = select_random_samples_with_weights(weights.clone(), 2);
        assert!(samples.is_ok());
        let samples = samples.unwrap();

        assert_eq!(samples.len(), 2);
        for &idx in &samples {
            assert!(idx < weights.len());
        }
    }

    #[test]
    fn test_full_sampling() {
        let weights = vec![0.5, 1.5, 2.0];
        let samples = select_random_samples_with_weights(weights.clone(), 3);
        assert!(samples.is_ok());
        let samples = samples.unwrap();

        // Must return all indices
        assert_eq!(samples, HashSet::from([0, 1, 2]));
    }

    #[test]
    fn test_heavy_weight_bias() {
        let weights = vec![1000.0, 0.0001, 0.0001];
        let mut counts = [0; 3];

        for _ in 0..100 {
            let s = select_random_samples_with_weights(weights.clone(), 1);
            assert!(s.is_ok());
            let s = s.unwrap();
            let idx = *s.iter().next().unwrap();
            counts[idx] += 1;
        }

        // Index 0 should dominate
        assert!(counts[0] > 90);
    }
}
