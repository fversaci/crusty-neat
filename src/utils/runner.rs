use super::config::RunConfiguration;
use super::fasta_tools::{read_fasta, write_fasta};
use super::fastq_tools::write_fastq;
use super::make_reads::generate_reads;
use super::mutate::mutate_fasta;
use super::nucleotides::Nuc;
use super::read_models::read_quality_score_model_json;
use super::vcf_tools::write_vcf;
use crate::utils::file_tools::check_create_dir;
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
    // check that the configuration is valid and create the output directory
    config.check()?;
    check_create_dir(config.output_dir.as_ref().unwrap())?;

    // Create the prefix of the files to write
    let output_prefix = format!(
        "{}/{}",
        config.output_dir.unwrap().display(),
        config.output_prefix.unwrap()
    );

    // Read the reference file into memory
    info!(
        "Mapping reference fasta file: {}",
        &config.reference.as_ref().unwrap().display()
    );
    let (fasta_map, fasta_order) = read_fasta(config.reference.as_ref().unwrap())?;

    // Load models that will be used for the runs.
    // For now we will use the one supplied, pulled directly from
    // NEAT2.0's original model.
    let default_quality_score_model_file = "models/neat_quality_score_model.json";
    let quality_score_model = read_quality_score_model_json(default_quality_score_model_file)?;

    // Mutate the reference and recording the variant locations.
    info!("Mutating reference.");
    let (mutated_map, variant_locations) = mutate_fasta(
        &fasta_map,
        config.minimum_mutations,
        config.mutation_rate,
        rng,
    )?;

    if config.produce_fasta == Some(true) {
        info!("Outputting fasta file");
        write_fasta(
            &mutated_map,
            &fasta_order,
            config.overwrite_output.unwrap(),
            &output_prefix,
        )?;
    }

    if config.produce_vcf == Some(true) {
        info!("Writing vcf file");
        write_vcf(
            &variant_locations,
            &fasta_order,
            config.ploidy.unwrap(),
            &config.reference.unwrap(),
            config.overwrite_output.unwrap(),
            &output_prefix,
            rng,
        )?;
    }

    if config.produce_fastq == Some(true) {
        let mut read_sets: HashSet<Vec<Nuc>> = HashSet::new();
        for (_name, sequence) in mutated_map.iter() {
            // defined as a set of read sequences that should cover
            // the mutated sequence `coverage` number of times
            let data_set = generate_reads(
                sequence,
                &config.read_len.unwrap(),
                &config.coverage.unwrap(),
                config.paired_ended.unwrap(),
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
            &output_prefix,
            config.overwrite_output.unwrap(),
            config.paired_ended.unwrap(),
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
    use super::*;
    use crate::create_rng;
    use tempdir::TempDir;

    #[test]
    fn test_runner() -> Result<()> {
        let tmp_dir = TempDir::new("crusty_neat")?;
        let config = RunConfiguration {
            reference: Some("test_data/H1N1.fa".into()),
            output_dir: Some(tmp_dir.path().to_path_buf()),
            output_prefix: Some("crusty_out".to_string()),
            overwrite_output: Some(false),
            produce_fasta: Some(true),
            ..Default::default()
        };
        let mut rng = create_rng(Some("Hello Cruel World"));
        run_neat(config, &mut rng)?;
        Ok(())
    }
}
