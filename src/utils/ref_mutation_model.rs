use crate::utils::mutation_model::{MutRateByContig, MutationModel, Region};
use crate::utils::types::SeqByContig;
use anyhow::{Result, anyhow};
use log::info;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

/// A model for genome mutations
#[derive(Serialize, Deserialize)]
pub struct RefMutationModel<'a> {
    /// Reference genome
    #[serde(skip)]
    pub ref_genome: Option<&'a SeqByContig>,
    /// Mutation rates by contig.
    pub mut_rates: Option<MutRateByContig>,
    /// Models for different mutation types.    
    pub mm: MutationModel,
}

impl<'a> RefMutationModel<'a> {
    /// set the same mutation rate for all, entire contigs
    fn set_default_mut_rates(&mut self, mutation_rate: Option<f64>) -> Result<()> {
        if mutation_rate.is_none() {
            return Err(anyhow!(
                "Mutation rate must be specified to set default mutation rates."
            ));
        }
        let mut mut_rates = MutRateByContig::new();
        for entry in self.ref_genome.unwrap().iter() {
            let (contig, seq) = (entry.key(), entry.value());
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
    pub fn new_all_contigs(genome: &'a SeqByContig, mutation_rate: f64) -> Result<Self> {
        let mm = MutationModel::default();
        let mut rmm = RefMutationModel {
            ref_genome: Some(genome),
            mut_rates: None,
            mm,
        };
        rmm.set_default_mut_rates(Some(mutation_rate))?;
        Ok(rmm)
    }
    /// Attach a reference genome to the mutation model.
    pub fn attach_ref_genome(
        &mut self,
        genome: &'a SeqByContig,
        def_mut_rate: Option<f64>,
    ) -> Result<()> {
        self.ref_genome = Some(genome);
        if self.mut_rates.is_none() {
            self.set_default_mut_rates(def_mut_rate)?;
        } else {
            self.check_ref_vs_mut_rates()?;
        }
        Ok(())
    }
    /// Check if reference genome is compatible with the mutation rates.
    pub fn check_ref_vs_mut_rates(&self) -> Result<()> {
        if self.mut_rates.is_none() {
            return Err(anyhow!("Mutation rates are not defined"));
        }
        let ref_genome = self.ref_genome.unwrap();
        for (contig, regions) in self.mut_rates.as_ref().unwrap() {
            if !ref_genome.contains_key(contig) {
                return Err(anyhow!("Contig {} not present in reference genome", contig));
            }
            for region in regions {
                if region.end > ref_genome.get(contig).unwrap().len() {
                    return Err(anyhow!("Region extends beyond contig {} length", contig));
                }
            }
        }
        Ok(())
    }
    /// Serialize mutation model to a yaml file.
    pub fn write_to_file(&self, output_prefix: &Path) -> Result<()> {
        let filename = output_prefix.with_extension("mut_model.yml");
        info!("Writing full mutation model to {}", filename.display());
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
