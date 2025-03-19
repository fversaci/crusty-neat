use anyhow::{Result, anyhow};
use flate2::{Compression, read::GzDecoder, write::GzEncoder};
use log::info;
use rand::Rng;
use rand::distr::Distribution;
use rand::distr::weighted::WeightedIndex;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read};
use std::path::Path;

/// A model for quality scores
#[derive(Serialize, Deserialize)]
pub struct QualityModel {
    /// first quality score
    zeroth: WeightedIndex<f64>,
    /// next quality score, based on read position and previous score
    next: Vec<Vec<WeightedIndex<f64>>>,
    /// max position assumed by the model, larger position values will
    /// be capped at max_pos
    #[serde(skip)]
    max_pos: usize,
}

/// probability densities (non-cumulative)
#[derive(Serialize, Deserialize)]
struct InputData {
    seeds: Vec<f64>,
    weights: Vec<Vec<Vec<f64>>>,
}

macro_rules! qs_char {
    ($q:expr) => {
        ($q + 33) as u8 as char
    };
}

impl QualityModel {
    /// Check if model params are ok and saves max position
    pub fn finalize(&mut self) -> Result<()> {
        let l1 = self.zeroth.weights().count();
        let l2 = self.next.first().unwrap().len();
        if l1 != l2 {
            return Err(anyhow!("Quality score model has incompatible dimensions"));
        }
        let no_prob = self
            .next
            .iter()
            .all(|v| v.iter().all(|wi| wi.weights().count() == l1));
        if no_prob {
            self.max_pos = l1;
            Ok(())
        } else {
            Err(anyhow!("Quality score model has incompatible dimensions"))
        }
    }
    /// Serialize mutation model to a gzipped yaml file.
    pub fn write_to_file(&self, output_prefix: &Path) -> Result<()> {
        let filename = output_prefix.with_extension("qs_model.yml.gz");
        info!("Writing quality score to {}", filename.display());
        let file = std::fs::File::create(&filename)?;
        let encoder = GzEncoder::new(file, Compression::default());
        let mut writer = BufWriter::new(encoder);
        serde_yaml::to_writer(&mut writer, self)?;
        Ok(())
    }
    /// Read mutation model from a yaml file
    pub fn from_file(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        let reader: Box<dyn Read> = if path.extension().is_some_and(|ext| ext == "gz") {
            Box::new(GzDecoder::new(file))
        } else {
            Box::new(file)
        };
        let reader = BufReader::new(reader);
        let mut conf: Self = serde_yaml::from_reader(reader)?;
        conf.finalize()?;
        Ok(conf)
    }
    /// read probability densities (non-cumulative) from json
    pub fn from_json(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let data: InputData = serde_json::from_reader(reader)?;

        let first = WeightedIndex::new(data.seeds).expect("Invalid probabilities");

        let next = data
            .weights
            .into_iter()
            .map(|v| {
                v.into_iter()
                    .map(|w| WeightedIndex::new(w).expect("Invalid probabilities"))
                    .collect()
            })
            .collect();
        let mut qs = Self {
            zeroth: first,
            next,
            max_pos: 0,
        };
        qs.finalize()?;
        Ok(qs)
    }
    /// get quality score for position 0
    fn get_first<R: Rng>(&self, rng: &mut R) -> usize {
        self.zeroth.sample(rng)
    }
    /// get next quality score
    fn get_next<R: Rng>(&self, pos: usize, prev: usize, rng: &mut R) -> usize {
        let wi = &self.next[pos][prev];
        wi.sample(rng)
    }
    pub fn generate_quality_scores<R: Rng>(
        &self,
        read_length: usize,
        rng: &mut R,
    ) -> Result<String> {
        let mut qs = String::with_capacity(read_length);
        // Insert the 0-th value
        let mut q = self.get_first(rng);
        qs.push(qs_char!(q));

        let waypoint = self.max_pos.min(read_length);

        // Generate quality scores up to min(max_pos, read_length)
        for p in 1..waypoint {
            q = self.get_next(p, q, rng);
            qs.push(qs_char!(q));
        }

        // Generate quality scores beyond max_pos (if needed)
        for _ in waypoint..read_length {
            q = self.get_next(waypoint, q, rng);
            qs.push(qs_char!(q));
        }

        Ok(qs)
    }
}
