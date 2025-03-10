use super::{
    mutation::MutationType,
    nucleotides::{self, Nuc},
};
use anyhow::{anyhow, Result};
use rand::distr::weighted::WeightedIndex;
use rand::distr::Distribution;
use rand::Rng;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Mutations independent of the reference genome.
#[derive(Serialize, Deserialize)]
pub struct MutationModel {
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
    pub fn default() -> Self {
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

impl InsModel {
    /// Get a mutated sequence
    pub fn get_alt_bases<R: Rng>(&self, ref_base: Nuc, rng: &mut R) -> Vec<Nuc> {
        let len = self.length.sample(rng);
        let mut seq = vec![ref_base];
        seq.extend(nucleotides::random_seq(rng, len));
        seq
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

impl DelModel {
    /// Get length of sequence to remove
    pub fn get_len<R: Rng>(&self, max_len: usize, rng: &mut R) -> usize {
        let len = self.length.sample(rng);
        len.min(max_len)
    }
}
