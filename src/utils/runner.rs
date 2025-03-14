use crate::utils::config::RunConfiguration;
use crate::utils::fasta_tools::{read_fasta, write_fasta};
use crate::utils::fastq_tools::write_fastq;
use crate::utils::file_tools::check_create_dir;
use crate::utils::make_reads::generate_reads;
use crate::utils::mutate::{apply_mutations, mutate_genome};
use crate::utils::nucleotides::Nuc;
use crate::utils::read_models::read_quality_score_model_json;
use crate::utils::ref_mutation_model::RefMutationModel;
use crate::utils::vcf_tools::write_vcf;
use anyhow::Result;
use log::{info, trace};
use rand::Rng;
use rand::seq::SliceRandom;
use std::collections::HashSet;

/// Mutate the reference genome and generate the output files
pub fn run_neat<R: Rng>(config: RunConfiguration, rng: &mut R) -> Result<()> {
    // check that the configuration is valid and create the output directory
    config.check()?;
    check_create_dir(config.output_dir.as_ref().unwrap())?;

    // Create the prefix of the files to write
    let output_prefix = config
        .output_dir
        .unwrap()
        .join(config.output_prefix.unwrap());

    // Read the reference file into memory
    info!(
        "Mapping reference fasta file: {}",
        config.reference.as_ref().unwrap().display()
    );
    let (ref_genome, contig_order) = read_fasta(config.reference.as_ref().unwrap())?;

    // Read mutation model from file or generate a default one
    let mut mut_model;
    if let Some(mutation_model_file) = &config.mutation_model {
        info!(
            "Reading mutation model from file {}",
            mutation_model_file.display()
        );
        mut_model = RefMutationModel::from_file(mutation_model_file)?;
        mut_model.attach_ref_genome(&ref_genome, config.def_mutation_rate)?;
    } else {
        info!("Creating mutation model from scratch");
        mut_model = RefMutationModel::new_all_contigs(
            &ref_genome,
            config.def_mutation_rate.unwrap_or(0.0),
        )?;
    }

    // serialize the mutation model to a yaml file
    info!(
        "Saving final mutation model to file {}",
        &output_prefix.display()
    );
    mut_model.write_to_file(&output_prefix)?;
    trace!("mutation model {:?}", mut_model.mm);

    // Mutate the reference and record the variant locations.
    info!("Mutating reference.");
    let mutations = mutate_genome(&ref_genome, &mut_model, rng)?;
    let mut_genome = apply_mutations(&ref_genome, &mutations)?;

    if config.produce_fasta == Some(true) {
        write_fasta(
            &mut_genome,
            &contig_order,
            config.overwrite_output.unwrap(),
            &output_prefix,
        )?;
    }

    if config.produce_vcf == Some(true) {
        write_vcf(
            &mutations,
            &contig_order,
            config.ploidy.unwrap(),
            mut_model.mm.prob_mut_multiple,
            &config.reference.unwrap(),
            config.overwrite_output.unwrap(),
            &output_prefix,
            rng,
        )?;
    }

    if config.produce_fastq == Some(true) {
        let mut read_sets: HashSet<Vec<Nuc>> = HashSet::new();
        for (_name, sequence) in mut_genome.iter() {
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

        // Load the quality score model, which simulates sequencing error rates
        // by generating quality scores for the produced reads.
        // Currently, we use the original model from NEAT2.0.
        let default_quality_score_model_file = "models/neat_quality_score_model.json";
        let quality_score_model = read_quality_score_model_json(default_quality_score_model_file)?;

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
    }
    info!("Processing complete");
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
