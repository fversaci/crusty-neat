// This is the run configuration for this particular run, which holds the parameters needed by the
// various side functions. It is build with a ConfigurationBuilder, which can take either a
// config yaml file or command line arguments and turn them into the configuration.

use super::cli::Cli;
use super::file_tools::check_create_dir;
use anyhow::{anyhow, Result};
use log::{info, warn};
use serde_yaml::Value;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::string::String;
use std::{env, fs};

#[derive(Debug)]
pub struct RunConfiguration {
    // This struct holds all the parameters for this particular run. It is derived from input either
    // from a configuration file or from command line inputs. This is not built directly in the code,
    // but is constructed by a builder to enable default values.
    //
    // reference: The path to the reference for the run.
    // read_len: The length of reads in the output fastq.
    // coverage: The average depth of coverage for the output fastq file.
    // mutation_rate: The rate of mutation for the file.
    // ploidy: The number of copies of each chromosome in the target organism. Mutation process will
    // be replicated this number of times.
    // paired_ended: If the run will be in paired-ended mode.
    // fragment_mean: Mean size of the fragments in paired-ended mode.
    // fragment_st_dev: Standard deviation of the fragment mean describing the sample set to sample
    // fragments from.
    // produce_fastq: True or false on whether to produce an output fastq file.
    // produce_fasta: True or false on whether to produce an output fasta file, 1 per ploid.
    // produce_vcf: True or false on whether to produce an output VCF file, with genotyped variants.
    // produce_bam: True or false on whether to produce an output BAM file, which will be aligned to
    // the reference.
    // overwrite_output: if true, will overwrite output. If false will error and exit you attempt to
    // overwrite files with the same name.
    // output_dir: The directory, relative or absolute, path to the directory to place output.
    // output_prefix: The name to use for the output files.
    pub reference: String,
    pub read_len: usize,
    pub coverage: usize,
    pub mutation_rate: f64,
    pub ploidy: usize,
    pub paired_ended: bool,
    pub fragment_mean: Option<f64>,
    pub fragment_st_dev: Option<f64>,
    pub produce_fastq: bool,
    pub produce_fasta: bool,
    pub produce_vcf: bool,
    pub produce_bam: bool,
    pub rng_seed: Option<String>,
    pub overwrite_output: bool,
    pub minimum_mutations: Option<usize>,
    pub output_dir: PathBuf,
    pub output_prefix: String,
}
#[allow(dead_code)]
impl RunConfiguration {
    // The purpose of this function is to redirect you to the ConfigBuilder
    pub fn build() -> Result<ConfigBuilder> {
        ConfigBuilder::new()
    }
}

// The config builder allows us to construct a config in multiple different ways, depending
// on the input.
pub struct ConfigBuilder {
    pub(crate) reference: Option<String>,
    read_len: usize,
    coverage: usize,
    mutation_rate: f64,
    ploidy: usize,
    paired_ended: bool,
    fragment_mean: Option<f64>,
    fragment_st_dev: Option<f64>,
    produce_fastq: bool,
    pub(crate) produce_fasta: bool,
    pub(crate) produce_vcf: bool,
    produce_bam: bool,
    rng_seed: Option<String>,
    overwrite_output: bool,
    pub(crate) minimum_mutations: Option<usize>,
    pub(crate) output_dir: PathBuf,
    output_prefix: String,
}

impl ConfigBuilder {
    pub fn new() -> Result<Self> {
        let cb = ConfigBuilder {
            // Setting default values
            reference: None,
            read_len: 150,
            coverage: 10,
            mutation_rate: 0.001,
            ploidy: 2,
            paired_ended: false,
            fragment_mean: None,
            fragment_st_dev: None,
            produce_fastq: true,
            produce_fasta: false,
            produce_vcf: false,
            produce_bam: false,
            rng_seed: None,
            overwrite_output: false,
            minimum_mutations: None,
            output_dir: env::current_dir()?,
            output_prefix: String::from("neat_out"),
        };
        Ok(cb)
    }

    pub fn check_and_print_config(&self) -> Result<()> {
        // This does a final check of the configuration for valid items. It will print info
        // message of the items, to work as a record and to assist in debugging any issues that
        // come up.

        let reference = self
            .reference
            .as_ref()
            .ok_or_else(|| anyhow!("No reference was specified."))?;
        info!(
            "Running rusty-neat to generate reads on {} with...",
            reference
        );

        info!("  >read length: {}", self.read_len);
        info!("  >coverage: {}", self.coverage);
        info!("  >mutation rate: {}", self.mutation_rate);
        info!("  >ploidy: {}", self.ploidy);
        info!("  >paired ended: {}", self.paired_ended);
        if self.overwrite_output {
            warn!("Overwriting any existing files.")
        }
        if let Some(min_mutations) = self.minimum_mutations {
            info!("  >minimum mutations per contig: {}", min_mutations);
        }
        let output_path = &self.output_dir;
        // This check may be overkill, but here it is. Let's make sure we ended up with something
        if !output_path.as_path().is_dir() {
            warn!(
                "Output directory is not a directory: {:?}",
                self.output_dir.display()
            );
            check_create_dir(output_path)?;
        }
        let file_prefix = format!("{}/{}", self.output_dir.display(), self.output_prefix);

        // No point in running if we aren't producing files
        if !(self.produce_fastq || self.produce_fasta || self.produce_vcf || self.produce_bam) {
            return Err(anyhow!(
                "All file types set to false, no files would be produced."
            ));
        }

        if self.paired_ended {
            if self.fragment_mean.is_none() || self.fragment_st_dev.is_none() {
                return Err(anyhow!(
                    "Paired ended is set to true, but fragment mean \
                    and standard deviation were not set."
                ));
            }
            if self.produce_fastq {
                info!("\t> fragment mean: {}", self.fragment_mean.unwrap());
                info!(
                    "\t> fragment standard deviation: {}",
                    self.fragment_st_dev.unwrap()
                );
                info!(
                    "Producing fastq files:\n\t> {}_r1.fastq\n\t {}_r2.fastq",
                    file_prefix, file_prefix
                )
            } else {
                info!("Producing fastq file:\n\t> {}_r1.fastq", file_prefix)
            }
        }
        if self.produce_fasta {
            info!("Producing fasta file: {}.fasta", file_prefix);
        }
        if self.produce_vcf {
            info!("Producing vcf file: {}.vcf", file_prefix)
        }
        if self.produce_bam {
            info!("Produce bam file: {}.bam", file_prefix)
        }
        if let Some(rng_seed) = self.rng_seed.as_ref() {
            info!("Using rng seed: {}", rng_seed);
        }
        Ok(())
    }

    // Function to build the actual configuration.
    pub fn build(self) -> Result<RunConfiguration> {
        let rc = RunConfiguration {
            reference: self.reference.ok_or_else(|| anyhow!("Missing reference"))?,
            read_len: self.read_len,
            coverage: self.coverage,
            mutation_rate: self.mutation_rate,
            ploidy: self.ploidy,
            paired_ended: self.paired_ended,
            fragment_mean: self.fragment_mean,
            fragment_st_dev: self.fragment_st_dev,
            produce_fastq: self.produce_fastq,
            produce_fasta: self.produce_fasta,
            produce_vcf: self.produce_vcf,
            produce_bam: self.produce_bam,
            rng_seed: self.rng_seed,
            overwrite_output: self.overwrite_output,
            minimum_mutations: self.minimum_mutations,
            output_dir: self.output_dir,
            output_prefix: self.output_prefix,
        };
        Ok(rc)
    }
}

fn generate_error(key: &str, key_type: &str, value: &Value) -> String {
    format!(
        "Input {} could not be converted to {}: {:?}",
        key, key_type, value
    )
}

pub fn read_config_yaml(yaml: String) -> Result<RunConfiguration> {
    // Reads an input configuration file from yaml using the serde package. Then sets the parameters
    // based on the inputs. A "." value means to use the default value.

    // Opens file for reading
    let file = fs::File::open(yaml)?;
    // Uses serde_yaml to read the file into a HashMap
    let scrape_config: HashMap<String, Value> = serde_yaml::from_reader(file)?;
    // Create the config builder then update any items from the configuration file in the
    // configuration object and returns it.
    let mut config_builder = ConfigBuilder::new()?;
    for (key, value) in scrape_config {
        // Two extra checks needed are for reference. Everything else can be
        // easily skipped with a value of "."
        if key == "reference" {
            if let Some(reference_str) = value.as_str() {
                let reference_path = Path::new(reference_str);
                if !reference_path.is_file() {
                    return Err(anyhow!("Reference file not found: {}", reference_str));
                } else {
                    config_builder.reference = reference_str.to_string().into();
                }
            } else {
                return Err(anyhow!(
                    "Expected a string for 'reference', but got: {:?}",
                    value
                ));
            }
        } else if value != "." {
            match key.as_str() {
                "read_len" => {
                    config_builder.read_len = value
                        .as_u64()
                        .ok_or_else(|| anyhow!(generate_error(&key, "integer", &value)))?
                        as usize
                }
                "coverage" => {
                    config_builder.coverage = value
                        .as_u64()
                        .ok_or_else(|| anyhow!(generate_error(&key, "integer", &value)))?
                        as usize
                }
                "mutation_rate" => {
                    config_builder.mutation_rate = value
                        .as_f64()
                        .ok_or_else(|| anyhow!(generate_error(&key, "float", &value)))?
                }
                "ploidy" => {
                    config_builder.ploidy = value
                        .as_u64()
                        .ok_or_else(|| anyhow!(generate_error(&key, "integer", &value)))?
                        as usize
                }
                "paired_ended" => {
                    config_builder.paired_ended = value
                        .as_bool()
                        .ok_or_else(|| anyhow!(generate_error(&key, "boolean", &value)))?
                }
                "fragment_mean" => {
                    config_builder.fragment_mean = value
                        .as_f64()
                        .ok_or_else(|| anyhow!(generate_error(&key, "float", &value)))?
                        .into() // to make it an option
                }
                "fragment_st_dev" => {
                    config_builder.fragment_st_dev = value
                        .as_f64()
                        .ok_or_else(|| anyhow!("{}", generate_error(&key, "float", &value)))?
                        .into() // to make it an option
                }
                "produce_fastq" => {
                    config_builder.produce_fastq = value
                        .as_bool()
                        .ok_or_else(|| anyhow!(generate_error(&key, "boolean", &value)))?
                }
                "produce_fasta" => {
                    config_builder.produce_fasta = value
                        .as_bool()
                        .ok_or_else(|| anyhow!(generate_error(&key, "boolean", &value)))?
                }
                "produce_vcf" => {
                    config_builder.produce_vcf = value
                        .as_bool()
                        .ok_or_else(|| anyhow!(generate_error(&key, "boolean", &value)))?
                }
                "produce_bam" => {
                    config_builder.produce_bam = value
                        .as_bool()
                        .ok_or_else(|| anyhow!(generate_error(&key, "boolean", &value)))?
                }
                "rng_seed" => {
                    config_builder.rng_seed = value.as_str().map(String::from);
                    // to make it an option
                }
                "overwrite_output" => {
                    config_builder.overwrite_output = value
                        .as_bool()
                        .ok_or_else(|| anyhow!(generate_error(&key, "boolean", &value)))?
                }
                "minimum_mutations" => {
                    config_builder.minimum_mutations = Some(
                        value
                            .as_u64()
                            .ok_or_else(|| anyhow!(generate_error(&key, "Valid integer", &value)))?
                            as usize,
                    )
                }
                "output_dir" => {
                    let output_path = value
                        .as_str()
                        .ok_or_else(|| anyhow!("Path must be a string"))?;
                    config_builder.output_dir = PathBuf::from(output_path);
                }
                "output_prefix" => {
                    config_builder.output_prefix = value
                        .as_str()
                        .ok_or_else(|| anyhow!("Prefix must be a string"))?
                        .to_string();
                }
                _ => continue,
            }
        }
    }
    let _ = &config_builder.check_and_print_config();
    config_builder.build()
}

pub fn build_config_from_args(args: Cli) -> Result<RunConfiguration> {
    // Takes in a bunch of args from a clap CLI and builds a config based on that. More CLI options
    // will need additional items entered here. To add them to the config, so they can be implemented.

    // Create the ConfigBuilder object with default values
    let mut config_builder = ConfigBuilder::new()?;
    // Can't do a run without a reference
    if !args.reference.is_empty() {
        config_builder.reference = args.reference.into();
    } else {
        return Err(anyhow!("No reference was specified."));
    }
    // The default value works directly for the config builder and CLI handles the type checking
    config_builder.read_len = args.read_length;
    config_builder.coverage = args.coverage;
    // default is empty string, in which case the config builder controls the default
    if args.output_dir.is_empty() {
        config_builder.output_dir = env::current_dir()?
    } else {
        let output_path = Path::new(&args.output_dir);
        check_create_dir(output_path)?;
        config_builder.output_dir = PathBuf::from(output_path);
    };
    // If this is unset, sets the default value of "neat_out" by CLI
    config_builder.output_prefix = args.output_file_prefix;
    // To set a minimum mutation rate, such as for debugging, or for small datasets, use this
    if let Some(min_mutations) = args.minimum_mutations {
        let input_min_muts = min_mutations as usize;
        config_builder.minimum_mutations = Some(input_min_muts);
    }
    let _ = &config_builder.check_and_print_config();
    config_builder.build()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_run_configuration() {
        let test_configuration = RunConfiguration {
            reference: String::from("Hello.world"),
            read_len: 100,
            coverage: 22,
            mutation_rate: 0.09,
            ploidy: 3,
            paired_ended: true,
            fragment_mean: Option::from(333.0),
            fragment_st_dev: Option::from(33.0),
            produce_fastq: false,
            produce_bam: true,
            produce_fasta: true,
            produce_vcf: true,
            rng_seed: None,
            overwrite_output: true,
            minimum_mutations: None,
            output_dir: PathBuf::from("/my/my"),
            output_prefix: String::from("Hey.hey"),
        };

        println!("{:?}", test_configuration);
        assert_eq!(test_configuration.reference, "Hello.world".to_string());
        assert_eq!(test_configuration.read_len, 100);
        assert_eq!(test_configuration.coverage, 22);
        assert_eq!(test_configuration.mutation_rate, 0.09);
        assert_eq!(test_configuration.ploidy, 3);
        assert!(test_configuration.paired_ended);
        assert_eq!(test_configuration.fragment_mean.unwrap(), 333.0);
        assert_eq!(test_configuration.fragment_st_dev.unwrap(), 33.0);
        assert!(!test_configuration.produce_fastq);
        assert!(test_configuration.produce_vcf);
        assert!(test_configuration.produce_bam);
        assert!(test_configuration.produce_fasta);
        assert_eq!(test_configuration.rng_seed, None);
        assert!(test_configuration.overwrite_output);
        assert_eq!(test_configuration.output_dir, PathBuf::from("/my/my"));
        assert_eq!(test_configuration.output_prefix, "Hey.hey".to_string());
    }

    #[test]
    fn test_build() -> Result<()> {
        use super::*;
        let x = RunConfiguration::build()?;
        assert!(x.reference.is_none());
        Ok(())
    }

    #[test]
    fn test_read_config_yaml() -> Result<()> {
        let yaml = String::from("config/neat_test.yml");
        let test_config = read_config_yaml(yaml)?;
        assert_eq!(test_config.reference, "test_data/ecoli.fa".to_string());
        assert_eq!(test_config.coverage, 3);
        Ok(())
    }

    #[test]
    fn test_bad_yaml() -> Result<()> {
        let yaml = String::from("fake_file.yml");
        let er = read_config_yaml(yaml);
        assert!(er.is_err());
        Ok(())
    }

    #[test]
    fn test_missing_ref() -> Result<()> {
        let yaml = String::from("config/simple_template.yml");
        let er = read_config_yaml(yaml);
        assert!(er.is_err());
        Ok(())
    }

    #[test]
    fn test_creates_out_dir() -> Result<()> {
        let yaml = String::from("config/neat_test_bad.yml");
        read_config_yaml(yaml)?;
        assert!(Path::new("fake").is_dir());
        Ok(fs::remove_dir("fake")?)
    }

    #[test]
    fn test_command_line_inputs() -> Result<()> {
        let args: Cli = Cli {
            config: String::new(),
            reference: String::from("test_data/ecoli.fa"),
            output_dir: String::from("test_data"),
            log_level: String::from("Trace"),
            log_dest: String::new(),
            output_file_prefix: String::from("test"),
            minimum_mutations: None,
            read_length: 150,
            coverage: 10,
        };

        let test_config = build_config_from_args(args)?;
        assert_eq!(test_config.reference, "test_data/ecoli.fa".to_string());
        Ok(())
    }

    #[test]
    fn test_bad_config_builder() -> Result<()> {
        let config = ConfigBuilder::new()?;
        let er = config.check_and_print_config();
        assert!(er.is_err());
        Ok(())
    }

    #[test]
    fn test_creat_nonexisting_out() -> Result<()> {
        let mut config = ConfigBuilder::new()?;
        config.reference = Some("test_data/H1N1.fa".to_string());
        config.output_dir = PathBuf::from("contig/");
        config.check_and_print_config()?;
        fs::remove_dir("contig")?;
        Ok(())
    }

    #[test]
    fn test_cl_missing_ref() -> Result<()> {
        let args: Cli = Cli {
            config: String::new(),
            reference: String::from(""),
            output_dir: String::from("test_dir"),
            log_level: String::from("Trace"),
            log_dest: String::new(),
            output_file_prefix: String::from("test"),
            minimum_mutations: None,
            read_length: 150,
            coverage: 10,
        };

        let er = build_config_from_args(args);
        assert!(er.is_err());
        Ok(())
    }

    #[test]
    fn test_overwrite_warn() -> Result<()> {
        let mut config = ConfigBuilder::new()?;
        config.reference = Option::from("test_data/H1N1.fa".to_string());
        config.overwrite_output = true;
        config.check_and_print_config()?;
        Ok(())
    }

    #[test]
    fn test_produce_fastq_messages() -> Result<()> {
        let mut config = ConfigBuilder::new()?;
        config.reference = Option::from("test_data/H1N1.fa".to_string());
        config.paired_ended = true;
        config.fragment_mean = Some(100.0);
        config.fragment_st_dev = Some(10.0);
        // tests first branch of if statement for paired_ended & produce_fastq = true
        config.check_and_print_config()?;
        // Checks the alternative pe = true, produce_fastq = false
        config.produce_fastq = false;
        // need to produce at least one file or check will raise an error
        config.produce_fasta = true;
        config.check_and_print_config()?;
        Ok(())
    }

    #[test]
    fn test_produce_fasta_messages() -> Result<()> {
        let mut config = ConfigBuilder::new()?;
        config.reference = Some("test_data/H1N1.fa".to_string());
        config.produce_fasta = true;
        config.check_and_print_config()?;
        config.produce_vcf = true;
        config.check_and_print_config()?;
        config.produce_bam = true;
        config.check_and_print_config()?;
        // If it passes all the checks, we're good.
        Ok(())
    }

    #[test]
    fn test_no_files() -> Result<()> {
        let mut config = ConfigBuilder::new()?;
        config.reference = Some("test_data/H1N1.fa".to_string());
        config.produce_fastq = false;
        let er = config.check_and_print_config();
        assert!(er.is_err());
        Ok(())
    }

    #[test]
    fn test_no_frag_mean_or_stdev() -> Result<()> {
        let mut config = ConfigBuilder::new()?;
        config.reference = Some("test_data/H1N1.fa".to_string());
        // paired end set to true, by default, fragment mean and st dev are None
        config.paired_ended = true;
        let er = config.check_and_print_config();
        assert!(er.is_err());
        Ok(())
    }

    #[test]
    fn test_no_frag_mean() -> Result<()> {
        let mut config = ConfigBuilder::new()?;
        config.reference = Some("test_data/H1N1.fa".to_string());
        config.paired_ended = true;
        config.fragment_st_dev = Some(10.0);
        let er = config.check_and_print_config();
        assert!(er.is_err());
        Ok(())
    }

    #[test]
    fn test_no_stdev() -> Result<()> {
        let mut config = ConfigBuilder::new()?;
        config.reference = Some("test_data/H1N1.fa".to_string());
        config.paired_ended = true;
        config.fragment_mean = Some(100.0);
        let er = config.check_and_print_config();
        assert!(er.is_err());
        Ok(())
    }

    #[test]
    fn no_output_dir_given() -> Result<()> {
        let args: Cli = Cli {
            config: String::new(),
            reference: String::from("test_data/H1N1.fa"),
            output_dir: String::new(),
            log_level: String::from("Trace"),
            log_dest: String::new(),
            output_file_prefix: String::from("test"),
            minimum_mutations: None,
            read_length: 150,
            coverage: 10,
        };

        let config = build_config_from_args(args)?;
        assert_eq!(env::current_dir()?.as_path(), config.output_dir);
        Ok(())
    }

    #[test]
    fn test_minimum_mutations_and_others() -> Result<()> {
        let args: Cli = Cli {
            config: String::new(),
            reference: String::from("test_data/H1N1.fa"),
            output_dir: String::new(),
            log_level: String::from("Trace"),
            log_dest: String::new(),
            output_file_prefix: String::from("test"),
            minimum_mutations: Some(10),
            read_length: 120,
            coverage: 13,
        };

        let config = build_config_from_args(args)?;
        assert_eq!(Some(10), config.minimum_mutations);
        assert_eq!(120, config.read_len);
        assert_eq!(13, config.coverage);
        Ok(())
    }
}
