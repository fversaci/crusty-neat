use super::quality_scores::QualityScoreModel;
use anyhow::Result;
use serde_json;
use std::fs;

/// Read a quality score model from a JSON file.
///
/// # Arguments
///
/// * `filename` - The name of the file to read the quality score model from.
pub fn read_quality_score_model_json(filename: &str) -> Result<QualityScoreModel> {
    let file = fs::File::open(filename)?;
    Ok(serde_json::from_reader(file)?)
}
