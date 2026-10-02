use std::collections::{HashMap, HashSet};
use serde::Deserialize;
use csv::ReaderBuilder;
use nori::load_factor_graph_bytes;

/// Represents a single effect weight record parsed from a CSV file.
#[derive(Deserialize)]
pub struct EffectWeight {
    pub sequence: String,
    pub score: f32,
    pub effect: usize,
}


/// Parses a CSV string into a vector of `EffectWeight` structs.
///
/// # Arguments
/// * `sequence_scores_csv` - A string containing CSV data for effect weights. The CSV
///   must include headers: `id, sequence, score, psms, effect, weight, log_weight`.
///
/// # Returns
/// Returns a `Result` containing a vector of `EffectWeight` structs if parsing succeeds.
///
/// # Errors
/// Returns an error if the CSV cannot be read, or if any record fails deserialization.
pub fn parse_effect_weights_csv(sequence_scores_csv: String) -> Result<Vec<EffectWeight>, Box<dyn std::error::Error>> {
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .from_reader(sequence_scores_csv.as_bytes());

    let mut sequence_scores = Vec::new();
    for record in rdr.deserialize() {
        let row: EffectWeight = record?;
        sequence_scores.push(row);
    }

    Ok(sequence_scores)
}


/// Generates a GraphML representation of a factor graph from a CSV string of effect weights.
///
/// # Arguments
/// * `sequence_scores_csv` - A string containing CSV data for effect weights.
///
/// # Returns
/// Returns a `Result` containing a GraphML string representation of the factor graph.
///
/// # Errors
/// Returns an error if CSV parsing fails or if any error occurs during graph construction.
pub fn generate_graph(sequence_scores_csv: String) -> Result<Vec<u8>, Box<dyn std::error::Error>> {

    let sequence_scores = parse_effect_weights_csv(sequence_scores_csv)?;

    let peptide_effect_graph = effect_weights_to_graphml(&sequence_scores);
    let factor_graph = load_factor_graph_bytes(&peptide_effect_graph)?;
    Ok(factor_graph)
}

pub fn effect_weights_to_graphml(effect_weights: &Vec<EffectWeight>) -> String {

    let mut xml = String::new();
    xml.push_str(r#"<?xml version="1.0" encoding="UTF-8"?>"#);
    xml.push_str(r#"<graphml xmlns="http://graphml.graphdrawing.org/xmlns">"#);
    xml.push_str(r#"  <key id="type" for="node"/>"#);
    xml.push_str(r#"  <key id="belief" for="node"/>"#);
    xml.push_str(r#"  <graph edgedefault="undirected">"#);
    
    // Keep unique sequence nodes
    let mut seen_sequences = HashSet::new();
    // Store first score seen for each sequence
    let mut sequence_scores: HashMap<&str, f32> = HashMap::new();

    for item in effect_weights {
        sequence_scores
            .entry(item.sequence.as_str())
            .or_insert(item.score);
    }

    // Sequence nodes
    for (sequence, score) in &sequence_scores {
        if seen_sequences.insert(*sequence) {
            xml.push_str(&format!(
                r#"    <node id="{}">
      <data key="type">input</data>
      <data key="belief">[{}, {}]</data>
    </node>
"#,
                sequence,
                1.0 - score,
                score
            ));
        }
    }

    // Effects nodes
    let mut seen_effects = HashSet::new();
    for item in effect_weights {
        if seen_effects.insert(item.effect) {
            xml.push_str(&format!(
                r#"    <node id="{}">
      <data key="type">output</data>
    </node>
"#,
                item.effect
            ));
        }
    }

    // Edges
    for item in effect_weights {
        xml.push_str(&format!(
            r#"    <edge source="{}" target="{}"/>
"#,
            item.sequence,
            item.effect
        ));
    }
    xml.push_str("  </graph>\n");
    xml.push_str("</graphml>\n");

    xml
}


#[cfg(test)]
mod tests {
    use super::*;

    fn sample_csv() -> String {
        "id,sequence,score,psms,effect,weight,log_weight
1,PEPTIDE1,0.8,3,100,0.5,-0.3
2,PEPTIDE2,0.6,3,100,0.4,-0.5
3,PEPTIDE3,0.9,3,200,0.7,-0.1"
            .to_string()
    }

    #[test]
    fn test_parse_effect_weights_csv() {
        let csv = sample_csv();
        let effects = parse_effect_weights_csv(csv).unwrap();
        assert_eq!(effects.len(), 3);
        assert!((effects[1].score - 0.6).abs() < 1e-6);
    }
}
