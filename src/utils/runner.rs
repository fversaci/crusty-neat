use crate::utils::nucleotides::Nuc;

use super::config::RunConfiguration;
use super::fasta_tools::{read_fasta, write_fasta};
use super::fastq_tools::write_fastq;
use super::make_reads::generate_reads;
use super::mutate::mutate_fasta;
use super::read_models::read_quality_score_model_json;
use super::vcf_tools::write_vcf;
use anyhow::Result;
use log::info;
use rand::seq::SliceRandom;
use rand::Rng;
use std::collections::HashSet;

/// The main function that runs the NEAT simulation.
/// This function will mutate the reference genome, and then produce
/// the output files
///
/// # Arguments
///
/// * `config` - The configuration for the run
/// * `rng` - The random number generator to use
///
/// # Returns
///
/// * `Result<()>` - A result that will be Ok(()) if the run was
///   successful
pub fn run_neat<R: Rng>(config: RunConfiguration, rng: &mut R) -> Result<()> {
    // Create the prefix of the files to write
    let output_file = format!("{}/{}", config.output_dir.display(), config.output_prefix);

    // Reading the reference file into memory
    info!("Mapping reference fasta file: {}", &config.reference);
    let (fasta_map, fasta_order) = read_fasta(&config.reference)?;

    // Load models that will be used for the runs.
    // For now we will use the one supplied, pulled directly from NEAT2.0's original model.
    let default_quality_score_model_file = "models/neat_quality_score_model.json";
    let quality_score_model = read_quality_score_model_json(default_quality_score_model_file)?;

    // Mutating the reference and recording the variant locations.
    info!("Mutating reference.");
    let (mutated_map, variant_locations) = mutate_fasta(&fasta_map, config.minimum_mutations, rng)?;

    if config.produce_fasta {
        info!("Outputting fasta file");
        write_fasta(
            &mutated_map,
            &fasta_order,
            config.overwrite_output,
            &output_file,
        )?;
    }

    if config.produce_vcf {
        info!("Writing vcf file");
        write_vcf(
            &variant_locations,
            &fasta_order,
            config.ploidy,
            &config.reference,
            config.overwrite_output,
            &output_file,
            rng,
        )?;
    }

    if config.produce_fastq {
        let mut read_sets: HashSet<Vec<Nuc>> = HashSet::new();
        for (_name, sequence) in mutated_map.iter() {
            // defined as a set of read sequences that should cover
            // the mutated sequence `coverage` number of times
            let data_set = generate_reads(
                sequence,
                &config.read_len,
                &config.coverage,
                config.paired_ended,
                config.fragment_mean,
                config.fragment_st_dev,
                rng,
            )?;

            read_sets.extend(data_set);
        }

        info!("Shuffling output fastq data");
        let outsets: Vec<&Vec<Nuc>> = read_sets.iter().collect();
        let mut outsets_order: Vec<usize> = (0..outsets.len()).collect();
        outsets_order.shuffle(rng);

        info!("Writing fastq");
        write_fastq(
            &output_file,
            config.overwrite_output,
            config.paired_ended,
            outsets,
            outsets_order,
            quality_score_model,
            rng,
        )?;
        info!("Processing complete")
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::config::ConfigBuilder;
    use super::*;
    use crate::create_rng;
    use std::fs;
    use std::path::PathBuf;

    #[test]
    fn test_runner() -> Result<()> {
        let mut config = RunConfiguration::build()?;
        config.reference = Some("test_data/H1N1.fa".to_string());
        // Because we are building this the wrong way, we need to manually create the output dir
        config.output_dir = PathBuf::from("test");
        fs::create_dir("test")?;
        let config = config.build()?;
        let mut rng = create_rng(Some("Hello Cruel World"));
        run_neat(config, &mut rng)?;
        fs::remove_dir_all("test")?;
        Ok(())
    }

    #[test]
    fn test_runner_files_messages() -> Result<()> {
        let mut config = ConfigBuilder::new()?;
        config.reference = Some("test_data/H1N1.fa".to_string());
        config.produce_fasta = true;
        config.produce_vcf = true;
        // Because we are building this the wrong way, we need to manually create the output dir
        config.output_dir = PathBuf::from("output");
        fs::create_dir("output")?;
        let config = config.build()?;
        let mut rng = create_rng(Some("Hello Cruel World"));
        run_neat(config, &mut rng)?;
        fs::remove_dir_all("output")?;
        Ok(())
    }
}
