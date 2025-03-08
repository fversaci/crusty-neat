use crate::LevelFilter;
use anyhow::{anyhow, Result};
use clap::Parser;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct Args {
    /// Path to the configuration file
    #[arg(short = 'c', long)]
    pub config_file: Option<PathBuf>,
    /// Log level: TRACE, DEBUG, INFO, WARN, ERROR
    #[arg(long, default_value_t = LevelFilter::Info)]
    pub log_level: LevelFilter,
    #[arg(long)]
    pub log_dest: Option<PathBuf>,
    #[clap(flatten)]
    pub config: Option<RunConfiguration>,
}

#[derive(Default, Clone, Debug, Serialize, Deserialize, Parser)]
pub struct RunConfiguration {
    /// Path to the reference genome
    #[arg(short = 'r', long)]
    pub reference: Option<PathBuf>,
    /// Length of the reads
    #[arg(long)]
    pub read_len: Option<usize>,
    /// Number of reads to generate
    #[arg(long)]
    pub coverage: Option<usize>,
    /// Default mutation rate
    #[arg(long)]
    pub def_mutation_rate: Option<f64>,
    /// Path to the mutation model
    #[arg(long)]
    pub mutation_model: Option<PathBuf>,
    /// Ploidy of the organism
    #[arg(long)]
    pub ploidy: Option<usize>,
    /// Generate paired ended reads
    #[arg(long)]
    pub paired_ended: Option<bool>,
    /// Mean of the fragment size
    #[arg(long)]
    pub fragment_mean: Option<f64>,
    /// Standard deviation of the fragment size
    #[arg(long)]
    pub fragment_st_dev: Option<f64>,
    /// Produce fastq files
    #[arg(long)]
    pub produce_fastq: Option<bool>,
    /// Produce fasta files
    #[arg(long)]
    pub produce_fasta: Option<bool>,
    /// Produce vcf files
    #[arg(long)]
    pub produce_vcf: Option<bool>,
    /// Produce bam files
    #[arg(long)]
    pub produce_bam: Option<bool>,
    /// Seed for the random number generator
    #[arg(long)]
    pub rng_seed: Option<String>,
    /// Overwrite output files
    #[arg(long)]
    pub overwrite_output: Option<bool>,
    #[arg(long, short = 'o')]
    /// Output directory
    pub output_dir: Option<PathBuf>,
    /// Prefix for the output files
    #[arg(long)]
    pub output_prefix: Option<String>,
}

impl RunConfiguration {
    /// Create a new configuration with some default values
    pub fn fill() -> Self {
        RunConfiguration {
            reference: None,
            read_len: Some(150),
            coverage: Some(10),
            def_mutation_rate: Some(0.001),
            mutation_model: None,
            ploidy: Some(2),
            paired_ended: Some(false),
            fragment_mean: None,
            fragment_st_dev: None,
            produce_fastq: Some(false),
            produce_fasta: Some(false),
            produce_vcf: Some(false),
            produce_bam: Some(false),
            rng_seed: None,
            overwrite_output: Some(false),
            output_dir: None,
            output_prefix: Some("crusty_out".to_string()),
        }
    }
    /// Override the values of the current configuration with the
    /// non-None values of another configuration
    pub fn override_with(&mut self, other: &RunConfiguration) -> Result<()> {
        // Serialize both structs into JSON `Value`
        let mut a_json: Value = serde_json::to_value(&mut *self)?;
        let b_json: Value = serde_json::to_value(other)?;

        // Iterate over the fields of `b` and override non-`None` values in `a`
        if let Some(b_fields) = b_json.as_object() {
            for (key, b_value) in b_fields {
                if b_value != &Value::Null {
                    // Only override if the value in `b` is not `null`
                    if let Some(a_field) = a_json.get_mut(key) {
                        *a_field = b_value.clone();
                    }
                }
            }
        }

        // Deserialize back into `RunConfiguration`
        *self = serde_json::from_value(a_json)?;
        Ok(())
    }
    /// Read a configuration from a yaml file
    pub fn from_file(path: &PathBuf) -> Result<Self> {
        let file = std::fs::File::open(path)?;
        let reader = std::io::BufReader::new(file);
        let conf: Self = serde_yaml::from_reader(reader)?;
        Ok(conf)
    }
    /// Check if the configuration for the run is valid
    pub fn check(&self) -> Result<()> {
        // mandatory fields
        if self.reference.is_none() {
            return Err(anyhow::anyhow!("Reference genome is required"));
        }
        if self.output_dir.is_none() {
            return Err(anyhow::anyhow!("Output directory is required"));
        }
        if self.overwrite_output.is_none() {
            return Err(anyhow!("Overwrite output should be set to true or false"));
        }
        // produce at least one output
        if self.produce_fastq.is_none()
            && self.produce_fasta.is_none()
            && self.produce_vcf.is_none()
            && self.produce_bam.is_none()
        {
            return Err(anyhow::anyhow!(
                "At least one output format must be selected"
            ));
        }
        // check that mutation rate is between 0 and 0.5
        if let Some(rate) = self.def_mutation_rate {
            if !(0.0..=0.5).contains(&rate) {
                return Err(anyhow!("Mutation rate should be between 0 and 0.5"));
            }
        }
        // paired ended reads require fragment size
        if self.paired_ended == Some(true)
            && (self.fragment_mean.is_none() || self.fragment_st_dev.is_none())
        {
            return Err(anyhow!(
                "Fragment mean and standard deviation are required for paired ended reads"
            ));
        }
        // ploidy is required for VCF output
        if self.produce_vcf == Some(true) && self.ploidy.is_none() {
            return Err(anyhow!("Ploidy is required for VCF output"));
        }
        // read length, coverage and paired ended are required for FASTQ output
        if self.produce_fastq == Some(true)
            && (self.read_len.is_none() || self.coverage.is_none() || self.paired_ended.is_none())
        {
            return Err(anyhow!(
                "Read length, coverage and paired ended are required for FASTQ output"
            ));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_override_with() {
        let mut a = RunConfiguration {
            reference: Some(PathBuf::from("ref.fa")),
            read_len: Some(150),
            coverage: Some(10),
            def_mutation_rate: Some(0.001),
            mutation_model: None,
            ploidy: Some(2),
            paired_ended: Some(false),
            fragment_mean: None,
            fragment_st_dev: None,
            produce_fastq: Some(false),
            produce_fasta: Some(false),
            produce_vcf: Some(false),
            produce_bam: Some(false),
            rng_seed: None,
            overwrite_output: Some(false),
            output_dir: None,
            output_prefix: Some("crusty_out".to_string()),
        };

        let b = RunConfiguration {
            reference: Some(PathBuf::from("ref2.fa")),
            read_len: None,
            coverage: Some(20),
            def_mutation_rate: None,
            mutation_model: Some(PathBuf::from("some_model.yml")),
            ploidy: None,
            paired_ended: Some(true),
            fragment_mean: Some(200.0),
            fragment_st_dev: Some(20.0),
            produce_fastq: None,
            produce_fasta: Some(true),
            produce_vcf: None,
            produce_bam: None,
            rng_seed: Some("1234".to_string()),
            overwrite_output: None,
            output_dir: Some(PathBuf::from("output")),
            output_prefix: None,
        };

        a.override_with(&b).unwrap();

        assert_eq!(a.reference, Some(PathBuf::from("ref2.fa")));
        assert_eq!(a.read_len, Some(150));
        assert_eq!(a.coverage, Some(20));
        assert_eq!(a.def_mutation_rate, Some(0.001));
        assert_eq!(a.mutation_model, Some(PathBuf::from("some_model.yml")));
        assert_eq!(a.ploidy, Some(2));
        assert_eq!(a.paired_ended, Some(true));
        assert_eq!(a.fragment_mean, Some(200.0));
        assert_eq!(a.fragment_st_dev, Some(20.0));
        assert_eq!(a.produce_fastq, Some(false));
        assert_eq!(a.produce_fasta, Some(true));
        assert_eq!(a.produce_vcf, Some(false));
        assert_eq!(a.produce_bam, Some(false));
        assert_eq!(a.rng_seed, Some("1234".to_string()));
        assert_eq!(a.overwrite_output, Some(false));
        assert_eq!(a.output_dir, Some(PathBuf::from("output")));
        assert_eq!(a.output_prefix, Some("crusty_out".to_string()));
    }

    #[test]
    fn test_check() -> Result<()> {
        let mut conf = RunConfiguration::default();
        assert!(conf.check().is_err());

        conf.reference = Some(PathBuf::from("ref.fa"));
        assert!(conf.check().is_err());

        conf.output_dir = Some(PathBuf::from("output"));
        assert!(conf.check().is_err());

        conf.produce_fastq = Some(false);
        conf.produce_fasta = Some(false);
        conf.produce_vcf = Some(false);
        conf.produce_bam = Some(false);
        assert!(conf.check().is_err());

        conf.produce_fasta = Some(true);
        assert!(conf.check().is_err());

        conf.overwrite_output = Some(false);
        assert!(conf.check().is_ok());

        conf.def_mutation_rate = Some(-0.2);
        assert!(conf.check().is_err());
        conf.def_mutation_rate = Some(0.6);
        assert!(conf.check().is_err());
        conf.def_mutation_rate = Some(0.001);
        assert!(conf.check().is_ok());
        conf.def_mutation_rate = None;
        assert!(conf.check().is_ok());

        conf.produce_fastq = Some(true);
        assert!(conf.check().is_err());

        conf.read_len = Some(150);
        conf.coverage = Some(10);
        conf.paired_ended = Some(false);
        assert!(conf.check().is_ok());

        conf.paired_ended = Some(true);
        assert!(conf.check().is_err());

        conf.fragment_mean = Some(200.0);
        assert!(conf.check().is_err());

        conf.fragment_st_dev = Some(200.0);
        assert!(conf.check().is_ok());

        conf.produce_vcf = Some(true);
        assert!(conf.check().is_err());

        conf.ploidy = Some(2);
        assert!(conf.check().is_ok());

        Ok(())
    }
}
