use crate::utils::{
    mutation::{Mutation, MutationType},
    nucleotides::{self, Nuc},
};
use crate::{idx_to_nuc, nuc_to_idx};
use anyhow::{Result, anyhow};
use rand::Rng;
use rand::distr::Distribution;
use rand::distr::weighted::WeightedIndex;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Mutations independent of the reference genome.
#[derive(Debug, Default, Serialize, Deserialize)]
pub struct MutationModel {
    /// Probabilities of different mutation types.
    pub mut_probabilities: MutProbabilities,
    /// WightedIndex to sample mutation types.
    #[serde(skip)]
    pub mut_type: Option<WeightedIndex<f64>>,
    /// SNP mutation probabilities.
    pub snp_model: SnpModel,
    /// Trinucleotide SNP mutation probabilities.
    pub tri_snp_model: Option<TriSnpModel>,
    /// Insertion mutation probabilities.
    pub ins_model: InsModel,
    /// Deletion mutation probabilities.
    pub del_model: DelModel,
    /// Probability that a mutation affect more homologous copies of
    /// the chromosome
    pub prob_mut_multiple: f64,
}

impl MutationModel {
    /// Get the probability that a given mutation type occurs.
    pub fn get_mut_probability(&self, mut_type: MutationType) -> f64 {
        *self.mut_probabilities.p.get(&mut_type).unwrap_or(&0.0)
    }
    /// init mut_type based on the mutation probabilities
    pub fn init_mut_type(&mut self) -> Result<()> {
        if self.mut_type.is_some() {
            return Ok(());
        }
        let weights: Vec<f64> = self.mut_probabilities.p.values().copied().collect();
        self.mut_type = Some(WeightedIndex::new(weights)?);
        Ok(())
    }
    /// get mutation type based on the mutation probabilities
    pub fn get_mut_type<R: Rng>(&self, rng: &mut R) -> Result<MutationType> {
        let mut_type_index = self
            .mut_type
            .as_ref()
            .ok_or_else(|| anyhow!("Run init_mut_type before get_mut_type"))?
            .sample(rng);
        match mut_type_index {
            0 => Ok(MutationType::Snp),
            1 => Ok(MutationType::Ins),
            2 => Ok(MutationType::Del),
            _ => unreachable!(),
        }
    }
    /// Create new snp mutation at given position
    pub fn create_snp<R: Rng>(&self, seq: &[Nuc], pos: usize, rng: &mut R) -> Result<Mutation> {
        let ref_base = seq[pos];
        let tri_snp_model = &self.tri_snp_model;
        let alt_base = if tri_snp_model.is_none() || pos == 0 || pos > seq.len() - 2 {
            let snp_model = &self.snp_model;
            snp_model.choose_new_nuc(ref_base, rng)?
        } else {
            // get also context
            let ref_bases = &seq[pos - 1..pos + 2];
            let tri_snp_model = tri_snp_model.as_ref().unwrap();
            tri_snp_model.choose_new_nuc(ref_bases, rng)?
        };
        // construct and return the mutation
        Mutation::new_snp(pos, ref_base, alt_base)
    }
    /// Create new insert mutation at given position
    pub fn create_ins<R: Rng>(&self, seq: &[Nuc], pos: usize, rng: &mut R) -> Result<Mutation> {
        let ins_model = &self.ins_model;
        let ref_base = seq[pos];
        let alt_bases = ins_model.get_alt_bases(ref_base, rng);
        // construct and return the mutation
        Mutation::new_ins(pos, ref_base, alt_bases)
    }
    /// Create new deletion mutation at given position
    pub fn create_del<R: Rng>(&self, seq: &[Nuc], pos: usize, rng: &mut R) -> Result<Mutation> {
        let del_model = &self.del_model;
        let max_len = seq.len() - pos - 1;
        let len = 1 + del_model.get_len(max_len, rng); // delete after current base
        let ref_seq = seq[pos..pos + len].to_vec();
        // construct and return the mutation
        Mutation::new_del(pos, ref_seq)
    }
    /// create new mutation
    pub fn create_mutation<R: Rng>(
        &self,
        m: MutationType,
        seq: &[Nuc],
        pos: usize,
        rng: &mut R,
    ) -> Result<Mutation> {
        match m {
            MutationType::Snp => self.create_snp(seq, pos, rng),
            MutationType::Ins => self.create_ins(seq, pos, rng),
            MutationType::Del => self.create_del(seq, pos, rng),
        }
    }
}

/// A type for storing mutation rates by contig.
pub type MutRateByContig = HashMap<String, Vec<Region>>;

/// A struct for storing a region of a contig with a mutation rate.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Region {
    pub start: usize,
    pub end: usize,
    pub rate: f64,
}

/// A struct for storing probabilities of different mutation types.
#[derive(Debug, Serialize, Deserialize)]
pub struct MutProbabilities {
    p: HashMap<MutationType, f64>,
}

impl Default for MutProbabilities {
    fn default() -> Self {
        let mut mp = MutProbabilities { p: HashMap::new() };
        mp.p.insert(MutationType::Snp, 0.5);
        mp.p.insert(MutationType::Ins, 0.25);
        mp.p.insert(MutationType::Del, 0.25);
        mp
    }
}

/// A model for SNP mutations.
#[derive(Debug, Serialize, Deserialize)]
pub struct SnpModel {
    weights: [WeightedIndex<f64>; 4],
}

impl Default for SnpModel {
    /// Default mutation model based on the original from NEAT 2.0
    fn default() -> Self {
        let w_a = WeightedIndex::new(vec![0.0, 0.17, 0.69, 0.14]).unwrap();
        let w_c = WeightedIndex::new(vec![0.16, 0.0, 0.17, 0.67]).unwrap();
        let w_g = WeightedIndex::new(vec![0.67, 0.17, 0.0, 0.16]).unwrap();
        let w_t = WeightedIndex::new(vec![0.14, 0.69, 0.17, 0.0]).unwrap();
        let weights = [w_a, w_c, w_g, w_t];
        Self { weights }
    }
}

impl SnpModel {
    /// Given a base, choose a new base based on the weights in the model
    pub fn choose_new_nuc<R: Rng>(&self, base: Nuc, rng: &mut R) -> Result<Nuc> {
        // Pick the weights list for the base that was input
        let dist = &self.weights[nuc_to_idx!(base)?];
        // Sample the distribution
        idx_to_nuc!(dist.sample(rng))
    }
}

/// A model for trinucleotide SNP mutations.
#[derive(Debug, Serialize, Deserialize)]
pub struct TriSnpModel {
    weights: [[[WeightedIndex<f64>; 4]; 4]; 4],
}

impl Default for TriSnpModel {
    /// Default mutation model based on the original from NEAT 2.0
    fn default() -> Self {
        let w_a = WeightedIndex::new(vec![0.0, 0.17, 0.69, 0.14]).unwrap();
        let w_c = WeightedIndex::new(vec![0.16, 0.0, 0.17, 0.67]).unwrap();
        let w_g = WeightedIndex::new(vec![0.67, 0.17, 0.0, 0.16]).unwrap();
        let w_t = WeightedIndex::new(vec![0.14, 0.69, 0.17, 0.0]).unwrap();
        let w_1d = [w_a, w_c, w_g, w_t];
        let w_2d = [w_1d.clone(), w_1d.clone(), w_1d.clone(), w_1d];
        let weights = [w_2d.clone(), w_2d.clone(), w_2d.clone(), w_2d];
        Self { weights }
    }
}

impl TriSnpModel {
    /// Given a trinucleotide, choose a new base based on the weights in the model
    pub fn choose_new_nuc<R: Rng>(&self, bases: &[Nuc], rng: &mut R) -> Result<Nuc> {
        if bases.len() != 3 {
            return Err(anyhow!("Expected three bases for trinucleotide mutation"));
        }
        let prev = nuc_to_idx!(bases[0])?;
        let base = nuc_to_idx!(bases[1])?;
        let next = nuc_to_idx!(bases[2])?;
        let wei = &self.weights[prev][next];
        let dist = &wei[base];
        // Sample the distribution
        idx_to_nuc!(dist.sample(rng))
    }
}

// A model for insertion mutations.
#[derive(Debug, Serialize, Deserialize)]
pub struct InsModel {
    /// Weighted index for the length of the insertion.
    length: WeightedIndex<f64>,
}

impl Default for InsModel {
    /// Default insertion model
    fn default() -> Self {
        Self {
            length: WeightedIndex::new(vec![0.0, 0.25, 0.25, 0.25, 0.25]).unwrap(),
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
#[derive(Debug, Serialize, Deserialize)]
pub struct DelModel {
    /// Weighted index for the length of the deletion.
    length: WeightedIndex<f64>,
}

impl Default for DelModel {
    /// Default deletion model
    fn default() -> Self {
        Self {
            length: WeightedIndex::new(vec![0.0, 0.25, 0.25, 0.25, 0.25]).unwrap(),
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
