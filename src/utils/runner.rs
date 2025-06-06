use std::collections::HashMap;

use crate::utils::config::RunConfiguration;
use crate::utils::fasta_tools::{read_fasta, write_fasta};
use crate::utils::fastq_tools::write_fastq;
use crate::utils::file_tools::check_create_dir;
use crate::utils::make_reads::generate_reads;
use crate::utils::mutate::{apply_mutations, mutate_genome};
use crate::utils::nucleotides::Nuc;
use crate::utils::quality_model::QualityModel;
use crate::utils::ref_mutation_model::RefMutationModel;
use crate::utils::vcf_tools::{VcfParams, write_vcf};
use anyhow::Result;
use log::{info, trace};
use rand::Rng;
use rand::seq::SliceRandom;

/// Mutate the reference genome and generate the output files
pub fn run_neat<R: Rng + Clone + Send + Sync>(config: RunConfiguration, rng: &mut R) -> Result<()> {
    // check that the configuration is valid and create the output directory
    config.check()?;
    check_create_dir(config.output_dir.as_ref().unwrap())?;

    // Create the prefix of the files to write
    let output_prefix = config
        .output_dir
        .unwrap()
        .join(config.output_prefix.unwrap());

    // Read the reference file into memory
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

    // compute HashMap of contig lengths
    let mut contig_lengths = HashMap::new();
    for entry in mut_genome.iter() {
        let contig = entry.key();
        let seq = entry.value();
        contig_lengths.insert(contig.to_string(), seq.len());
    }

    if config.produce_vcf == Some(true) {
        let vcf_params = VcfParams {
            contig_lengths,
            mutations: &mutations,
            contig_order: &contig_order,
            ploidy: config.ploidy.unwrap(),
            prob_mut_multiple: mut_model.mm.prob_mut_multiple,
            reference_path: &config.reference.unwrap(),
            overwrite_output: config.overwrite_output.unwrap(),
            output_prefix: &output_prefix,
        };
        write_vcf(vcf_params, rng)?;
    }

    if config.produce_fastq == Some(true) {
        info!("Generating reads");
        // Convert DashMap to Vec<Vec<Nuc>>
        let sequences: Vec<Vec<Nuc>> = mut_genome.into_iter().map(|(_, v)| v).collect();
        let mut reads: Vec<&[Nuc]> = Vec::new();
        for sequence in &sequences {
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

            reads.extend(data_set);
        }

        info!("Shuffling reads");
        reads.shuffle(rng);

        // Load the quality score model, which simulates sequencing
        // error rates by generating quality scores for the produced
        // reads. Currently, we use the original model from NEAT2.0.
        // Read mutation model from file or generate a default one
        let quality_model;
        if let Some(quality_model_file) = &config.quality_model {
            info!(
                "Reading quality model from file {}",
                quality_model_file.display()
            );
            quality_model = QualityModel::from_file(quality_model_file)?;
        } else {
            return Err(anyhow::anyhow!("Quality model file not provided"));
        };

        info!("Writing reads to fastq");
        write_fastq(
            &output_prefix,
            config.overwrite_output.unwrap(),
            config.paired_ended.unwrap(),
            reads,
            &quality_model,
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
