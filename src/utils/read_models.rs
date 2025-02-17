use super::quality_scores::QualityScoreModel;
use anyhow::Result;
use serde_json;
use std::fs;

pub fn read_quality_score_model_json(filename: &str) -> Result<QualityScoreModel> {
    let file = fs::File::open(filename)?;
    Ok(serde_json::from_reader(file)?)
}
