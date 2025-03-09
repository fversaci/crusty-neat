use super::{
    mutation::Mutation,
    mutation_model::{MutationModel, Region},
    nucleotides::Nuc,
    types::{MutByContig, SeqByContig},
};
use anyhow::{anyhow, Result};
use log::debug;
use rand::Rng;
use rand_distr::{Distribution, Binomial};
use std::collections::{HashMap, HashSet};

/// Mutates a genome by generating mutations for each contig.
///
/// # Arguments
///
/// * `genome` - A hashmap of contig names (keys) and a vector
///   representing the original sequence.
/// * `mutation_model` - The mutation model for the run
/// * `rng` - random number generator for the run
///
/// # Returns
///
/// A tuple with pointers to:
/// * A vector of mutations for each contig
pub fn mutate_genome<'a, R: Rng>(
    genome: &'a SeqByContig,
    mutation_model: &MutationModel,
    rng: &mut R,
) -> Result<MutByContig<'a>> {
    let mut mutations: MutByContig = HashMap::new();
    let mut_rates = mutation_model.mut_rates.clone().unwrap();
    for (contig, regions) in mut_rates {
        let sequence = &genome[&contig];
        for region in regions {
            debug!(
                "Segment of {}: from {} to {}.",
                contig, region.start, region.end
            );
            // Generate mutations for the region
            let contig_mutations = mutate_region(sequence, mutation_model, &region, rng)?;
            // Add the generated mutations to the list of mutations for this contig
            mutations
                .entry(contig.clone())
                .or_default()
                .extend(contig_mutations);
        }
    }
    Ok(mutations)
}

/// This function takes a slice of Nuc's and mutates num_mutations at
/// random. It returns the vector of mutations.
///
/// # Arguments
///
/// * `sequence` - A vector representing the original sequence.
/// * `num_mutations` - The number of mutations to add to this
///   sequence.
/// * `rng` - random number generator for the run
///
/// # Returns
///
/// A tuple with pointers to:
/// * A vector of mutations for this region
fn mutate_region<'a, R: Rng>(
    sequence: &'a [Nuc],
    mutation_model: &MutationModel,
    region: &Region,
    rng: &mut R,
) -> Result<Vec<Mutation<'a>>> {
    // Calculate the number of mutations using a Binomial distribution
    let region_len = region.end - region.start;
    let bin = Binomial::new(region_len as u64, region.rate)?;
    let num_mutations: usize = bin.sample(rng).try_into()?;
    debug!("Adding {} mutations", num_mutations);
    let mut indexes_to_mutate: HashSet<usize> = HashSet::new();
    // choose num_positions distinct indexes to mutate
    while indexes_to_mutate.len() < num_mutations {
        let index = rng.random_range(region.start..region.end);
        if sequence[index] == Nuc::N {
            continue;
        }
        indexes_to_mutate.insert(index);
    }
    // Get the snp model
    let snp_model = &mutation_model.snp_model;
    // Will hold the variants added to this sequence
    let mut sequence_variants: Vec<Mutation<'a>> = Vec::new();
    // for each index, picks a new base
    for pos in indexes_to_mutate {
        // remember the reference for later.
        let reference_base = sequence[pos];
        // pick a new base and assign the position to it.
        let alt_base = snp_model.choose_new_nuc(reference_base, rng)?;
        // construct the snp mutation
        let snp = Mutation::new_snp(sequence, pos, alt_base)?;
        sequence_variants.push(snp);
    }

    Ok(sequence_variants)
}

/// Applies the mutations to the reference genome and returns the new genome.
pub fn apply_mutations<'a>(
    genome: &'a SeqByContig,
    mutations: &MutByContig<'a>,
) -> Result<SeqByContig> {
    let mut new_genome = genome.clone();
    for (contig, mutations) in mutations {
        let sequence = new_genome.get_mut(contig).unwrap();
        for mutation in mutations {
            match mutation {
                Mutation::Snp { pos, alt_base, .. } => {
                    sequence[*pos] = *alt_base;
                }
                _ => return Err(anyhow!("Non-SNP mutations not supported yet")),
            }
        }
    }
    Ok(new_genome)
}

#[cfg(test)]
mod tests {
    use anyhow::anyhow;

    use super::*;
    use crate::{create_rng, utils::nucleotides::random_seq};

    #[test]
    fn test_mutate_genome() -> Result<()> {
        let mut rng = create_rng(Some("Hello Cruel World"));
        let seq: Vec<Nuc> = random_seq(&mut rng, 100);
        let genome: SeqByContig = HashMap::from([("chr1".to_string(), seq.clone())]);
        let mut_model = MutationModel::all_contigs(&genome, 0.1)?;
        let mutations = mutate_genome(&genome, &mut_model, &mut rng)?;
        assert!(mutations.contains_key("chr1"));
        let mutation = mutations["chr1"][0].clone();
        match mutation {
            Mutation::Snp {
                pos,
                ref_base,
                alt_base,
            } => {
                assert_eq!(seq[pos], *ref_base);
                assert_ne!(alt_base, *ref_base);
            }
            _ => return Err(anyhow!("Expected a SNP mutation")),
        }
        Ok(())
    }

    #[test]
    fn test_mutate_genome_no_mutations() -> Result<()> {
        let mut rng = create_rng(Some("Hello Cruel World"));
        let seq: Vec<Nuc> = random_seq(&mut rng, 100);
        let genome: SeqByContig = HashMap::from([("chr1".to_string(), seq.clone())]);
        let mut_model = MutationModel::all_contigs(&genome, 0.0)?;
        let mutations = mutate_genome(&genome, &mut_model, &mut rng)?;
        assert!(mutations.contains_key("chr1"));
        assert!(mutations["chr1"].is_empty());
        Ok(())
    }
}
