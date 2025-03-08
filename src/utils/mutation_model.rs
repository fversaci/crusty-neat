use super::{mutation::MutationType, nucleotides::Nuc, types::SeqByContig};
use anyhow::{anyhow, Result};
use log::info;
use rand::distr::weighted::WeightedIndex;
use rand::distr::Distribution;
use rand::Rng;
use serde::{Deserialize, Serialize};
use std::{collections::HashMap, path::PathBuf};

/// A model for mutation rates and probabilities of different mutation types.
#[derive(Serialize, Deserialize)]
pub struct MutationModel {
    /// Mutation rates by contig.
    pub mut_rates: Option<MutRateByContig>,
    /// Probabilities of different mutation types.
    pub mut_probabilities: MutProbabilities,
    /// SNP mutation probabilities.
    pub snp_model: SnpModel,
    /// Insertion mutation probabilities.
    pub ins_model: InsModel,
    /// Deletion mutation probabilities.
    pub del_model: DelModel,
}

impl MutationModel {
    /// Get the probability that a given mutation type occurs.
    pub fn get_mut_probability(&self, mut_type: MutationType) -> f64 {
        *self.mut_probabilities.p.get(&mut_type).unwrap_or(&0.0)
    }
    pub fn set_default_mut_rates(&mut self, genome: &SeqByContig, mutation_rate: Option<f64>) -> Result<()> {
        if mutation_rate.is_none() {
            return Err(anyhow!("Mutation rate must be specified to set default mutation rates."));
        }
        let mut mut_rates = MutRateByContig::new();
        for (contig, seq) in genome {
            let mut_rates_by_contig = vec![Region {
                start: 0,
                end: seq.len(),
                rate: mutation_rate.unwrap(),
            }];
            mut_rates.insert(contig.to_string(), mut_rates_by_contig);
        }
        let mut_rates = Some(mut_rates);
        self.mut_rates = mut_rates;
        Ok(())
    }
    /// Create a mutation model with the given mutation rate for all contigs.
    pub fn all_contigs(genome: &SeqByContig, mutation_rate: f64) -> Result<Self> {
        let snp_model = SnpModel::default();
        let ins_model = InsModel::default();
        let del_model = DelModel::default();
        let mut mm = MutationModel {
            mut_rates: None,
            mut_probabilities: MutProbabilities::default(),
            snp_model,
            ins_model,
            del_model,
        };
        mm.set_default_mut_rates(genome, Some(mutation_rate))?;
        Ok(mm)
    }
    // Serialize mutation model to a yaml file.
    pub fn write_to_file(&self, output_prefix: &str) -> Result<()> {
        let filename = format!("{}.mut_model.yml", output_prefix);
        info!("Writing full mutation model to {}", filename);
        let file = std::fs::File::create(&filename)?;
        serde_yaml::to_writer(&file, self)?;
        Ok(())
    }
    /// Read mutation model from a yaml file
    pub fn from_file(path: &PathBuf) -> Result<Self> {
        let file = std::fs::File::open(path)?;
        let reader = std::io::BufReader::new(file);
        let conf: Self = serde_yaml::from_reader(reader)?;
        Ok(conf)
    }
}

/// A type for storing mutation rates by contig.
pub type MutRateByContig = HashMap<String, Vec<Region>>;

/// A struct for storing a region of a contig with a mutation rate.
#[derive(Clone, Serialize, Deserialize)]
pub struct Region {
    pub start: usize,
    pub end: usize,
    pub rate: f64,
}

/// A struct for storing probabilities of different mutation types.
#[derive(Serialize, Deserialize)]
pub struct MutProbabilities {
    p: HashMap<MutationType, f64>,
}

impl MutProbabilities {
    fn default() -> Self {
        let mut mp = MutProbabilities { p: HashMap::new() };
        mp.p.insert(MutationType::Snp, 0.5);
        mp.p.insert(MutationType::Ins, 0.25);
        mp.p.insert(MutationType::Del, 0.25);
        mp
    }
}

/// A model for SNP mutations.
#[derive(Serialize, Deserialize)]
pub struct SnpModel {
    a: WeightedIndex<f64>,
    c: WeightedIndex<f64>,
    g: WeightedIndex<f64>,
    t: WeightedIndex<f64>,
}

impl Default for SnpModel {
    /// Default mutation model based on the original from NEAT 2.0
    fn default() -> Self {
        Self {
            a: WeightedIndex::new(vec![0.0, 0.17, 0.69, 0.14]).unwrap(),
            c: WeightedIndex::new(vec![0.16, 0.0, 0.17, 0.67]).unwrap(),
            g: WeightedIndex::new(vec![0.67, 0.17, 0.0, 0.16]).unwrap(),
            t: WeightedIndex::new(vec![0.14, 0.69, 0.17, 0.0]).unwrap(),
        }
    }
}

impl SnpModel {
    /// Given a base, choose a new base based on the weights in the model
    pub fn choose_new_nuc<R: Rng>(&self, base: Nuc, rng: &mut R) -> Result<Nuc> {
        // Pick the weights list for the base that was input
        let dist = match base {
            Nuc::A => &self.a,
            Nuc::C => &self.c,
            Nuc::G => &self.g,
            Nuc::T => &self.t,
            _ => return Err(anyhow!("Invalid input base")),
        };
        // Sample the distribution
        match dist.sample(rng) {
            0 => Ok(Nuc::A),
            1 => Ok(Nuc::C),
            2 => Ok(Nuc::G),
            3 => Ok(Nuc::T),
            _ => Err(anyhow!("Invalid output base")),
        }
    }
}

/// A model for insertion mutations.
#[derive(Serialize, Deserialize)]
pub struct InsModel {
    /// Weighted index for the length of the insertion.
    length: WeightedIndex<f64>,
}

impl Default for InsModel {
    /// Default insertion model
    fn default() -> Self {
        Self {
            length: WeightedIndex::new(vec![0.25, 0.25, 0.25, 0.25]).unwrap(),
        }
    }
}

/// A model for deletion mutations.
#[derive(Serialize, Deserialize)]
pub struct DelModel {
    /// Weighted index for the length of the deletion.
    length: WeightedIndex<f64>,
}

impl Default for DelModel {
    /// Default deletion model
    fn default() -> Self {
        Self {
            length: WeightedIndex::new(vec![0.25, 0.25, 0.25, 0.25]).unwrap(),
        }
    }
}
